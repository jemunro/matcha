# matcha

## Project Description

matcha is a compiled, efficient structural variant (SV) matching and annotation tool written in Nim using hts-nim. It addresses limitations in existing tools (svafotate, truvari anno, SURVIVOR) which are slow, memory-intensive, or have awkward dependencies.

matcha has four planned modes:

- `matcha match` — find all pairwise matches between two SV callsets — **implemented (milestone 1)**
- `matcha anno` — annotate a query callset against a population database VCF, transferring arbitrary INFO fields — *planned*
- `matcha collapse` — merge SVs from multiple callers within a single sample into a consensus callset — *planned*
- `matcha merge` — merge SVs across multiple samples into a cohort-level pVCF — *planned*

The core matching engine is shared across all modes. The primary use case is short-read SV callsets (Manta, DELLY etc), so sequence-based matching is out of scope. All matching is coordinate and size based.

## Status

**Milestone 1 (`matcha match`)** — complete.

- Source modules in [src/matcha/](src/matcha/) all green.
- 44 tests passing across [tests/test_intervals.nim](tests/test_intervals.nim) (15), [tests/test_preproc.nim](tests/test_preproc.nim) (14), [tests/test_match.nim](tests/test_match.nim) (15).
- Inputs accepted as `.vcf.gz` or `.bcf` (format auto-detected from magic bytes).
- Output is tab-separated with a single `#`-prefixed header line.
- Verbose logging via `-v` / `--verbose` (top-level or after the subcommand) writes timestamped progress to stderr.

Build & run: `nimble build` produces `./matcha`. Run all tests: `nimble test`. Regenerate fixtures: `python3 tests/generate_fixtures.py` (requires `bcftools`).

---

## matcha match — Overview

Finds all pairwise matches between records in callsetA and callsetB. Only matched pairs are emitted. Operates on DEL, DUP, INV SVTYPE only (BND and INS are skipped silently). Genotypes are ignored.

### Matching logic

- Two records are candidates if they share CHROM and SVTYPE, and their positions are within a dynamic window
- Window size is derived analytically from SVLEN and the minimum threshold: `window = SVLEN * (1 - threshold)`
- Two metrics are computed for each candidate pair:
  - **Reciprocal overlap**: `overlap / min(lenA, lenB)` where overlap = `min(END_A, END_B) - max(POS_A, POS_B)`
  - **Jaccard**: `overlap / union` where union = `max(END_A, END_B) - min(POS_A, POS_B)`
- A pair is emitted if both metrics meet their respective minimums (unspecified metric defaults to 0.0)
- At least one of `--min-overlap` or `--min-jaccard` must be supplied

### Output

Tab-separated to stdout or `--output` file. The first line is a `#`-prefixed header so downstream tools can skip it (`awk '!/^#/'`, `grep -v ^#`). Columns:

```
#CHROM  POS_A  END_A  ID_A  POS_B  END_B  ID_B  SVTYPE  OVERLAP  JACCARD
```

### CLI

```
matcha [-v|--verbose] match [options] callsetA callsetB

Inputs may be VCF.gz (.vcf.gz) or BCF (.bcf); format is auto-detected.

Options:
  --min-overlap FLOAT             minimum reciprocal overlap (0.0-1.0)
  --min-jaccard FLOAT             minimum Jaccard index (0.0-1.0)
  --threads INT                   number of worker threads (default: 1)
  --tmp-dir PATH                  temp directory (default: system temp)
  --output PATH                   output file (default: stdout)
  -v, --verbose                   verbose logging to stderr
  -h, --help                      show this help
```

- One of `--min-overlap` or `--min-jaccard` is required
- `-v` / `--verbose` is accepted at the top level (`matcha -v match …`) or after the subcommand (`matcha match -v …`)

---

## Architecture

### Internal preprocessing

On invocation, matcha:

1. Reads both inputs (VCF.gz or BCF), normalizes each record (see semantics below), and writes temporary per-SVTYPE BCFs. One BCF per supported SVTYPE per input callset (so up to 6 temp files for DEL/DUP/INV across two callsets), each spanning all chromosomes.
2. Indexes each temp BCF with CSI.
3. Tracks which chromosomes appear per SVTYPE so the work queue can be built without re-scanning.
4. Builds a work queue of `(chrom, svtype)` pairs present in both callsets. Both jobs of a given SVTYPE share the same per-SVTYPE BCF; the worker uses a region query to scope to its chromosome.

#### Preprocessing semantics

Each record is resolved, validated, and slimmed before it lands in the per-SVTYPE temp BCF:

| Field | Resolution rule | Failure |
|---|---|---|
| **SVTYPE** | Prefer `INFO/SVTYPE`. Fall back to symbolic `ALT` (`<DEL>`/`<DUP>`/`<INV>`/`<INS>`/`<TRA>`). On disagreement, **ALT wins**. | Unresolvable → warn + skip. Resolved but not in `{DEL, DUP, INV}` → silent count under `unsupported_svtype` (BND/INS/TRA awaiting later milestone). |
| **END** | Prefer `INFO/END`. Fall back to `POS + abs(INFO/SVLEN)`. | Neither → warn + skip (`unresolvable_end`). After resolution, `END ≤ POS` → warn + skip (`end_le_pos`). |
| **SVLEN** | Prefer `INFO/SVLEN` (stored as `abs`). Fall back to `END − POS`. | Neither → warn + skip. If both `INFO/END` and `INFO/SVLEN` were independently provided and they disagree by >10%, warn (no skip) and use `END − POS` as authoritative. |
| **ID** | Synthesize `CHROM_POS_SVTYPE_LINENUMBER` (1-based input order) when absent or `.`. | Always succeeds. |
| **Contig** | `BCF_ERR_CTG_UNDEF` from hts-nim → warn + skip (`unknown_contig`). | — |

The output BCF carries only the **keep-set** INFO fields: `SVTYPE`, `SVLEN`, `END`, `CHR2`, `END2`, `POS2`. All others (DP, AF, AC, FORMAT data, samples, etc.) are dropped. The header still defines other lines that came in via `copy_header` (slimming the writer header would require `bcf_translate`, which hts-nim does not bind), but record bodies are minimal.

Warnings are emitted to stderr with the prefix `[matcha preproc WARN]`. Per-record warnings are throttled at 5 per reason per callset (override with `MATCHA_WARN_CAP`); an end-of-callset summary always emits, listing read/kept/skipped counts, per-reason counters, inconsistencies, and synthesized IDs.

### Matching

A thread pool of size `--threads` pulls jobs from the work queue. Each job:

1. Streams records from the A BCF restricted to its chromosome (`vcfA.query(job.chrom)`)
2. For each record, computes the query window from SVLEN
3. Queries the B BCF index for candidates in `[POS - window, END + window]`
4. Computes reciprocal overlap and Jaccard for each candidate
5. Writes matches to a per-job temp TSV

The thread pool uses a single shared atomic counter; workers `fetchAdd` to claim the next job index. With `--threads 1` the pool is bypassed and jobs run inline on the main thread.

### Output assembly

Main thread concatenates per-job TSVs in a consistent order (sorted by CHROM then SVTYPE) to the output destination, then removes the temp TSVs and BCFs.

---

## Module Structure

```
matcha/
  matcha.nimble        # nimble package definition
  config.nims          # vendored hts-nim path; panics:off for tests
  src/
    matcha.nim         # binary entry; `include matcha/main`
    matcha/
      main.nim         # CLI parsing (parseopt), subcommand dispatch
      utils.nim        # shared types: SvType, MatchResult, MatchConfig
      intervals.nim    # reciprocalOverlap, jaccard, queryWindow
      preproc.nim      # VCF/BCF → per-SVTYPE temp BCF split, CSI indexing, work queue
      match.nim        # per-job matching, thread pool, output assembly
      log.nim          # verbose logging (stderr, timestamped)
  tests/
    fixtures/          # generated by tests/generate_fixtures.py
      expected/        # expected TSVs for default / strict / jaccard_only
    generate_fixtures.py
    test_utils.nim     # `timed` template + watchdog timeout
    test_intervals.nim
    test_preproc.nim
    test_match.nim
  vendor/
    hts-nim/           # vendored
    blocky/
```

---

## Testing

### Fixtures

Fixtures are produced by [tests/generate_fixtures.py](tests/generate_fixtures.py) (Python; requires `bcftools` and `bgzip` on PATH). It writes:

- `tests/fixtures/fixtureA.vcf.gz`, `fixtureB.vcf.gz` (+ `.csi` indexes)
- `tests/fixtures/fixtureA.bcf`, `fixtureB.bcf` (+ `.csi` indexes) — same content, BCF format, used to test BCF input support
- `tests/fixtures/expected/expected_default.tsv` — `--min-overlap 0.5`
- `tests/fixtures/expected/expected_strict.tsv` — `--min-overlap 0.8 --min-jaccard 0.8`
- `tests/fixtures/expected/expected_jaccard_only.tsv` — `--min-jaccard 0.5`

The fixtures cover: exact match, partial overlap (above and below threshold), no overlap, SVTYPE mismatch, multiple matches, unmatched A record, size asymmetry, multi-chromosome (chr1/chr2/chrX), all three SVTYPEs, and BND/INS records that should be silently skipped.

### Running tests

```
nimble test                                      # build + all three suites
nim c --hints:off -r tests/test_intervals.nim    # one suite
```

Each test prints `<ID>\tPASS\t<elapsed>\t<desc>` or `<ID>\tFAIL\t…`. Per-test timeout is enforced by the `timed` template ([tests/test_utils.nim](tests/test_utils.nim)) — 10s default, override with `MATCHA_TEST_TIMEOUT`.

### Coverage

- **test_intervals.nim** — pure metric math: exact / partial / no overlap, adjacency, degenerate intervals, size asymmetry, query-window edge cases.
- **test_preproc.nim** — per-SVTYPE BCF split, `chromsBySvtype` tracking, BND/INS exclusion, field extraction, CSI presence, work-queue intersection, `.bcf` input support, INFO slim to keep-set, ID synthesis, ALT-fallback / ALT-wins-on-conflict SVTYPE resolution, skip rules (missing END/SVLEN, END ≤ POS), END/SVLEN inconsistency normalization.
- **test_match.nim** — end-to-end matcha invocations against fixtures: each TC from the fixture set, threshold filtering, multi-chromosome output, output column shape and `#`-prefixed header, `--output` file mode, multi-thread parity with single-thread, `.bcf` vs `.vcf.gz` parity.

---

## Dependencies

- [`hts-nim`](vendor/hts-nim/) — VCF/BCF reading, index querying *(vendored)*
- `bcftools` (test fixture generation only)
- Standard library: `os`, `parseopt`, `strutils`, `tables`, `atomics`, `algorithm`, `times`

CLI parsing uses `std/parseopt`, not `cligen`.

---

## Future Milestones (out of scope for milestone 1)

- BND / TRA matching
- INS matching
- `matcha anno` mode
- `matcha collapse` mode
- `matcha merge` mode
- Genotype handling
- Sequence similarity
