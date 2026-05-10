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
- 60 tests passing across [tests/test_intervals.nim](tests/test_intervals.nim) (15), [tests/test_bins.nim](tests/test_bins.nim) (14), [tests/test_preproc.nim](tests/test_preproc.nim) (16), [tests/test_match.nim](tests/test_match.nim) (15).
- Inputs accepted as `.vcf.gz` or `.bcf` (format auto-detected from magic bytes).
- Output is tab-separated with a single `#`-prefixed header line.
- Verbose logging via `-v` / `--verbose` (top-level or after the subcommand) writes timestamped progress to stderr.

Build & run: `nimble build` produces `./matcha`. Run all tests: `nimble test`. Regenerate fixtures: `python3 tests/generate_fixtures.py` (requires `bcftools`).

---

## matcha match — Overview

Finds all pairwise matches between records in callsetA and callsetB. Only matched pairs are emitted. Operates on DEL, DUP, INV SVTYPE only (BND and INS are skipped silently). Genotypes are ignored.

### Matching logic

- Two records are candidates if they share CHROM and SVTYPE and their sizes are within a ratio bounded by the threshold
- Two metrics are computed for each candidate pair:
  - **Reciprocal overlap**: `overlap / max(lenA, lenB)` where overlap = `min(END_A, END_B) - max(POS_A, POS_B)`. Equivalent to requiring `overlap/lenA ≥ t AND overlap/lenB ≥ t` simultaneously (standard truvari/bedtools definition).
  - **Jaccard**: `overlap / union` where union = `max(END_A, END_B) - min(POS_A, POS_B)`
- A pair is emitted if both supplied metrics meet their respective minimums; an unset metric is not checked
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

1. Reads both inputs (VCF.gz or BCF), normalizes each record (see semantics below), assigns it a log2 size bin, and writes it to a temporary per-(SVTYPE, bin) BCF. Up to ~30 temp files per callset (one per populated (svtype, bin) pair), each spanning all chromosomes.
2. Indexes each temp BCF with CSI.
3. Tracks which chromosomes appear per SVTYPE and which bins are populated per SVTYPE, and captures chromosome order from the input header.
4. Builds a work queue of `(chrom, svtype, binA)` jobs. Each job covers one chromosome, one SVTYPE, and one A size-bin. The adjacent B bins are precomputed from the binding threshold (stricter of `--min-overlap`/`--min-jaccard`) and the canonical reciprocal-overlap size constraint; jobs with no adjacent populated B bins are silently skipped.

#### Preprocessing semantics

Each record is resolved, validated, and slimmed before it lands in the per-(svtype, bin) temp BCF:

| Field | Resolution rule | Failure |
|---|---|---|
| **SVTYPE** | Prefer `INFO/SVTYPE`. Fall back to symbolic `ALT` (`<DEL>`/`<DUP>`/`<INV>`/`<INS>`/`<TRA>`). On disagreement, **ALT wins**. | Unresolvable → warn + skip. Resolved but not in `{DEL, DUP, INV}` → silent count under `unsupported_svtype` (BND/INS/TRA awaiting later milestone). |
| **END** | Prefer `INFO/END`. Fall back to `POS + abs(INFO/SVLEN)`. | Neither → warn + skip (`unresolvable_end`). After resolution, `END ≤ POS` → warn + skip (`end_le_pos`). |
| **SVLEN** | Prefer `INFO/SVLEN` (stored as `abs`). Fall back to `END − POS`. | Neither → warn + skip. If both `INFO/END` and `INFO/SVLEN` were independently provided and they disagree by >10%, warn (no skip) and use `END − POS` as authoritative. |
| **ID** | Synthesize `CHROM_POS_SVTYPE_LINENUMBER` (1-based input order) when absent or `.`. | Always succeeds. |
| **Contig** | `BCF_ERR_CTG_UNDEF` from hts-nim → warn + skip (`unknown_contig`). | — |
| **MATCHA_BOFF** | *Added* (not preserved): the BGZF virtual offset of the source-file record (`(block_address << 16) \| block_offset`), encoded as `Number=2 Type=Integer` (high32, low32). matcha-internal — points back into the *original* input file. | Always succeeds. |
| **Size bin** | `binIndexFor(svlen)` (log2 from 1024bp; bin 0 = `[0,1024)`). Determines which per-(svtype, bin) BCF receives the record. | Clamped to 0 for non-positive SVLEN; cannot fail. |

The output BCF carries only the **keep-set** INFO fields: `SVTYPE`, `SVLEN`, `END`, `CHR2`, `END2`, `POS2`, `MATCHA_BOFF`. All others (DP, AF, AC, FORMAT data, samples, etc.) are dropped. The header still defines other lines that came in via `copy_header` (slimming the writer header would require `bcf_translate`, which hts-nim does not bind), but record bodies are minimal.

Temp BCF names follow the pattern `matcha_<pid>_<A|B>_<SVTYPE>_b<bin>.bcf`. Writers are opened lazily — only (svtype, bin) pairs that receive at least one record get a file.

Warnings are emitted to stderr with the prefix `[matcha preproc WARN]`. Per-record warnings are throttled at 5 per reason per callset (override with `MATCHA_WARN_CAP`); an end-of-callset summary always emits, listing read/kept/skipped counts, per-reason counters, inconsistencies, and synthesized IDs.

### Size bins

SV lengths are partitioned on a log2 scale:

- **Bin 0**: `[0, 1024)` bp — catch-all for small SVs; `SVLEN < 1024` (and the impossible `SVLEN ≤ 0` defence).
- **Bin N** (N ≥ 1): `[2^(N+9), 2^(N+10))` bp — so bin 1 = `[1024, 2048)`, bin 2 = `[2048, 4096)`, etc.

`binIndexFor(svlen)` uses `fastLog2` from `std/bitops`.

#### Adjacent bins

For canonical reciprocal overlap (`overlap/max(lenA, lenB) ≥ t`), a B record of length `L_b` can only match an A record of length `L_a` if `L_b / L_a ≥ t` and `L_a / L_b ≥ t`, i.e. `L_b ∈ [L_a·t, L_a/t]`. Given A is in bin with range `[L_lo, L_hi)`, the eligible B range is `(L_lo·t, L_hi/t)`. `adjacentBins(binA, threshold, populatedB)` returns the B bins from `populatedB` whose range intersects this interval.

The threshold used for adjacency is `max(minOverlap, minJaccard)` (the binding constraint, since both metrics share the same size-ratio requirement). The per-pair filter still tests each metric individually.

### Matching

A thread pool of size `--threads` pulls jobs from the work queue. Each job (`runMatchJob`):

1. Opens the per-(svtype, binA) A BCF and streams records restricted to `job.chrom`.
2. For each A record, queries each adjacent B bin through a **TiledBuffer** — a lazy, eviction-based cache keyed by tile index.
3. Computes reciprocal overlap and Jaccard for each candidate; emits `MatchResult` rows that pass both thresholds.
4. Returns a `seq[MatchResult]` carrying both A's and B's `MATCHA_BOFF` source-file offsets.

**Position window**: asymmetric `[posA − U, posA + svlenA)` where `U` = tile width (upper bound of the B bin's size range). A B record up to `U` bp left of `posA` can extend rightward into A; a B record at `posA + svlenA` or later cannot overlap.

**Tiled buffer**: each B bin is cached in fixed-width tiles of width `U` (the B bin's upper bound). Tiles are loaded lazily on first access via a CSI region query. Eviction: tile `K` is safe to drop once `posA ≥ (K+2)·U`, because the earliest future A' needs tile `K` only if `posA' < (K+2)·U` (derivation: `queryStart = posA' − U < (K+1)·U`). CSI region queries can return records straddling tile boundaries; the fetcher assigns each record to the tile containing its POS to avoid double-counting.

The thread pool uses a single shared atomic counter; workers `fetchAdd` to claim the next job index and write their result into a disjoint slot in `gJobResults: seq[seq[MatchResult]]`. With `--threads 1` the pool is bypassed and jobs run inline on the main thread.

When `--threads ≥ 2`, preprocessing of A and B runs in two parallel threads before the work queue is built.

### Output assembly

Main thread iterates the per-job result slots in sorted order (CHROM then SVTYPE) and writes formatted rows to the output destination, prefixed by the `#`-header line. Temp BCFs and CSI indexes are removed at the end of the run.

### Source-record retrieval (anno / merge / collapse)

Preprocessing records each surviving record's `bgzf_tell` virtual offset into `INFO/MATCHA_BOFF` on the slim BCF. Downstream modes that need to read the original record (e.g. `anno` transferring user INFO fields) collect offsets from the match results, sort them ascending, then `bgzf_seek` + `bcf_read` in order. Sorted seeks decompress each touched BGZF block at most once — same amortisation as a region query, and strictly better for sparse retrievals. `bgzf_seek` itself is one extra binding (not present in vendored hts-nim, will be added locally when the first consumer lands).

Same-source exclusion (needed for `collapse`/`merge` to suppress self-matches) is not in scope for milestone 1. When that mode lands, add an optional `excludeFn: proc(idA, idB: string): bool` field to `MatchConfig` (nil = no-op); `runMatchJob` checks it before appending to results.

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
      bins.nim         # log2 size-bin assignment, adjacentBins, TiledBuffer
      preproc.nim      # VCF/BCF → per-(SVTYPE, bin) temp BCF split, CSI indexing, work queue
      match.nim        # per-job matching (bins + tiled buffers), thread pool, output assembly
      log.nim          # verbose logging (stderr, timestamped)
  tests/
    fixtures/          # generated by tests/generate_fixtures.py
      expected/        # expected TSVs for default / strict / jaccard_only
    generate_fixtures.py
    test_utils.nim     # `timed` template + watchdog timeout
    test_intervals.nim
    test_bins.nim      # B-prefix: bin math; T-prefix: TiledBuffer
    test_preproc.nim
    test_match.nim
  vendor/
    hts-nim/           # git submodule, pinned to v0.3.31
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

The fixtures cover: exact match, partial overlap (above and below threshold), no overlap, SVTYPE mismatch, multiple matches, unmatched A record, size asymmetry (which now scores `overlap/max` so large-vs-small pairs score low), multi-chromosome (chr1/chr2/chrX), all three SVTYPEs, BND/INS records that should be silently skipped, and at least one large-SV pair (≥1024bp) that lands in a non-zero size bin.

### Running tests

```
nimble test                                      # build + all four suites
nim c --hints:off -r tests/test_intervals.nim    # one suite
```

Each test prints `<ID>\tPASS\t<elapsed>\t<desc>` or `<ID>\tFAIL\t…`. Per-test timeout is enforced by the `timed` template ([tests/test_utils.nim](tests/test_utils.nim)) — 10s default, override with `MATCHA_TEST_TIMEOUT`.

### Coverage

- **test_intervals.nim** — pure metric math: exact / partial / no overlap, adjacency, degenerate intervals, size asymmetry (`overlap/max` denominator), query-window edge cases.
- **test_bins.nim** (B-prefix) — `binIndexFor` for SVLENs spanning bin 0 through large bins; `binRange` round-trip; `adjacentBins` at t=0.5 and t≈1.0; exclusion of non-adjacent bins. (T-prefix) — `TiledBuffer`: cold fetch, re-query with no extra fetcher calls, empty-tile memoization, spanning-two-tiles A record, eviction preserving current and just-passed tiles.
- **test_preproc.nim** — per-(svtype, bin) BCF split, `populatedBins` tracking, `chromsBySvtype` tracking, BND/INS exclusion, field extraction, CSI presence, work-queue intersection by `(chrom, svtype, binA)`, VCF header chrom order preserved in job ordering, `.bcf` input support, INFO slim to keep-set, ID synthesis, ALT-fallback / ALT-wins-on-conflict SVTYPE resolution, skip rules (missing END/SVLEN, END ≤ POS), END/SVLEN inconsistency normalization, BGZF seek round-trip, large-SV record landing in non-zero bin.
- **test_match.nim** — end-to-end matcha invocations against fixtures: each TC from the fixture set, threshold filtering (including `DEL_A_08` dropping under canonical RO), multi-chromosome output, output column shape and `#`-prefixed header, `--output` file mode, multi-thread parity with single-thread, `.bcf` vs `.vcf.gz` parity.

---

## Dependencies

- [`hts-nim`](vendor/hts-nim/) — VCF/BCF reading, index querying *(git submodule under `vendor/`, pinned to v0.3.31)*
- `bcftools` (test fixture generation only)
- Standard library: `os`, `parseopt`, `strutils`, `tables`, `sets`, `atomics`, `algorithm`, `bitops`, `times`

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
