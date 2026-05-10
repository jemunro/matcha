# matcha

## Project Description

matcha is a compiled, efficient structural variant (SV) matching and annotation tool written in Nim using hts-nim. It has four planned modes:

- `matcha match` — find all pairwise matches between two SV callsets — **implemented (milestone 1)**
- `matcha anno` — annotate a query callset against a population database VCF, transferring arbitrary INFO fields — *planned*
- `matcha collapse` — merge SVs from multiple callers within a single sample into a consensus callset — *planned*
- `matcha merge` — merge SVs across multiple samples into a cohort-level pVCF — *planned*

The core matching engine is shared across all modes. All matching is coordinate and size based; sequence similarity and genotypes are out of scope.

## Status

**Milestone 1 (`matcha match`)** — complete.

Build: `nimble build` → `./matcha`. Test: `nimble test`. Regenerate fixtures: `python3 tests/generate_fixtures.py` (needs `bcftools`).

---

## matcha match — Overview

Finds all pairwise matches between records in callsetA and callsetB. Only matched pairs are emitted. Operates on DEL, DUP, INV SVTYPE only (BND and INS are skipped silently). Genotypes are ignored.

### Matching logic

- Two records are candidates if they share CHROM and SVTYPE and their sizes are within a ratio bounded by the threshold.
- Two metrics are computed for each candidate pair:
  - **Reciprocal overlap**: `overlap / max(lenA, lenB)` where overlap = `min(END_A, END_B) - max(POS_A, POS_B)`. Equivalent to requiring `overlap/lenA ≥ t AND overlap/lenB ≥ t` simultaneously (standard truvari/bedtools definition).
  - **Jaccard**: `overlap / union` where union = `max(END_A, END_B) - min(POS_A, POS_B)`.
- A pair is emitted if both supplied metrics meet their respective minimums; an unset metric is not checked.
- At least one of `--min-overlap` or `--min-jaccard` must be supplied.

### Self-matching mode (`--self`)

`matcha match --self INPUT` matches a single callset against itself. The engine reuses the cross-callset pipeline; the only differences are a single preprocessing pass (`filesB` aliases `filesA`) and a per-pair filter `aOff < bOff` on `MATCHA_BOFF`, which both deduplicates the symmetric `(X,Y)`/`(Y,X)` ordering and excludes the trivial self-self case. The work queue prunes `binsB` to `{b ≥ binA}` so each cross-bin pair is built once.

### Output

Tab-separated to stdout or `--output` file. The first line is a `#`-prefixed header so downstream tools can skip it (`awk '!/^#/'`, `grep -v ^#`). Columns:

```
#CHROM  POS_A  END_A  ID_A  POS_B  END_B  ID_B  SVTYPE  OVERLAP  JACCARD
```

### CLI

```
matcha [-v|--verbose] match [options] callsetA callsetB    # cross-callset
matcha [-v|--verbose] match --self [options] INPUT         # self-match

Inputs may be VCF.gz (.vcf.gz) or BCF (.bcf); format is auto-detected.

Options:
  --min-overlap FLOAT             minimum reciprocal overlap (0.0-1.0)
  --min-jaccard FLOAT             minimum Jaccard index (0.0-1.0)
  --self                          match a single input against itself
                                  (each pair emitted once; no self-self)
  --threads INT                   number of worker threads (default: 1)
  --tmp-dir PATH                  temp directory (default: system temp)
  --output PATH                   output file (default: stdout)
  -v, --verbose                   verbose logging to stderr
  -h, --help                      show this help
```

- One of `--min-overlap` or `--min-jaccard` is required.
- `-v` / `--verbose` is accepted at the top level (`matcha -v match …`) or after the subcommand.
- `--self` requires exactly one positional input; passing two is an error.

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

The output BCF carries only the **keep-set** INFO fields: `SVTYPE`, `SVLEN`, `END`, `CHR2`, `END2`, `POS2`, `MATCHA_BOFF`. All others (DP, AF, AC, FORMAT data, samples, etc.) are dropped.

Temp BCF names follow the pattern `matcha_<pid>_<A|B>_<SVTYPE>_b<bin>.bcf`. Writers are opened lazily — only (svtype, bin) pairs that receive at least one record get a file.

Warnings are emitted to stderr with the prefix `[matcha preproc WARN]`. Per-record warnings are throttled at 5 per reason per callset (override with `MATCHA_WARN_CAP`); an end-of-callset summary always emits.

### Size bins

SV lengths are partitioned on a log2 scale:

- **Bin 0**: `[0, 1024)` bp — catch-all for small SVs.
- **Bin N** (N ≥ 1): `[2^(N+9), 2^(N+10))` bp — so bin 1 = `[1024, 2048)`, bin 2 = `[2048, 4096)`, etc.

#### Adjacent bins

A and B can match only if their sizes are within ratio `t` (the binding threshold = `max(minOverlap, minJaccard)`). `adjacentBins(binA, threshold, populatedB)` returns the populated B bins whose size range satisfies that constraint; the per-pair filter still tests each metric individually. Derivation lives in [src/matcha/bins.nim](src/matcha/bins.nim).

### Matching

A thread pool of size `--threads` pulls jobs from the work queue. Each job opens its per-(svtype, binA) A BCF, streams records restricted to `job.chrom`, and queries each adjacent B bin through a **TiledBuffer** — a lazy, eviction-based cache keyed by tile index.

**Position window**: asymmetric `[posA − U, posA + svlenA)` where `U` = upper bound of the B bin's size range. **Tiled buffer**: each B bin is cached in fixed-width tiles of width `U`; loaded lazily via CSI region query, evicted once A advances past the point where any future A record could need them. CSI queries can return records straddling tile boundaries — the fetcher assigns each record to the tile containing its POS to avoid double-counting. Derivations live in [src/matcha/bins.nim](src/matcha/bins.nim).

The thread pool uses a single shared atomic counter; workers `fetchAdd` to claim the next job index and write their result into a disjoint slot. With `--threads 1` the pool is bypassed and jobs run inline. With `--threads ≥ 2`, preprocessing of A and B runs in two parallel threads before the work queue is built.

### Output assembly

Main thread iterates the per-job result slots in sorted order (CHROM then SVTYPE) and writes formatted rows to the output destination, prefixed by the `#`-header line. Temp BCFs and CSI indexes are removed at the end of the run.

### Source-record retrieval (anno / merge / collapse)

Each surviving record's `bgzf_tell` virtual offset is written to `INFO/MATCHA_BOFF` on the slim BCF. Downstream modes that need the original record collect offsets from match results, sort ascending, then `bgzf_seek` + `bcf_read` in order — same amortisation as a region query, strictly better for sparse retrievals. `bgzf_seek` will need a one-line local binding when the first consumer lands.

---

## Module structure

Sources under [src/matcha/](src/matcha/):

- `main.nim` — CLI parsing (`std/parseopt`), subcommand dispatch
- `utils.nim` — shared types: `SvType`, `MatchResult`, `MatchConfig`
- `intervals.nim` — `reciprocalOverlap`, `jaccard`
- `bins.nim` — log2 size-bin assignment, `adjacentBins`, `TiledBuffer`
- `preproc.nim` — VCF/BCF normalize → per-(svtype, bin) BCF + work queue
- `match.nim` — per-job matching, thread pool, output assembly
- `log.nim` — verbose logging (stderr, timestamped)

Tests under [tests/](tests/) mirror the module names; `test_utils.nim` provides the `timed` template + watchdog timeout. Vendored hts-nim lives in `vendor/hts-nim/` (submodule, pinned v0.3.31).

---

## Testing

Fixtures are generated by [tests/generate_fixtures.py](tests/generate_fixtures.py) (Python; needs `bcftools`+`bgzip` on PATH) into `tests/fixtures/`. Expected TSVs cover the default / strict / jaccard-only / self thresholds against the same fixtures.

Run all tests with `nimble test`. Each test prints `<ID>\tPASS|FAIL\t<elapsed>\t<desc>`; per-test timeout via the `timed` template (10s default, override with `MATCHA_TEST_TIMEOUT`).

---

## Dependencies

- [hts-nim](vendor/hts-nim/) — VCF/BCF I/O (git submodule, pinned v0.3.31)
- `bcftools` — fixture generation only
- Nim ≥ 2.0 standard library (CLI uses `std/parseopt`, not `cligen`).
