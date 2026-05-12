# matcha — Design & Architecture

Reference for contributors and Claude Code. For user documentation see [README.md](README.md).

---

## Pipeline overview

Three phases per invocation:

1. **Preproc** — normalize + slim each input into per-(SVTYPE, bin) temp BCFs, indexed with CSI.
2. **Match** — thread pool processes a work queue of `(chrom, svtype, binA)` jobs; each job pairs A records against adjacent B bins using TiledBuffer (intervals) or a sliding deque (BND).
3. **Output** — main thread collects per-job results in sorted order and writes TSV (`match`) or annotated VCF/BCF (`anno`).

---

## Key invariants

- **SVTYPE resolution**: ALT wins over `INFO/SVTYPE` on disagreement; ALT is authoritative for BND mate (`CHR2`/`POS2`).
- **MATCHA_BOFF**: BGZF virtual offset of the source record, stored as `Number=2 Integer` (high32, low32). Used by anno (and future modes) to `bgzf_seek` back to the original record.
- **Self-mode dedup**: `aOff < bOff` filter on `MATCHA_BOFF` eliminates symmetric pairs and self-self. Interval queue prunes `binsB` to `{b ≥ binA}`.
- **BND always in bin 0**: BND records are point events, not size-binned. They always land in `(svBND, 0)` temp BCFs and bypass `adjacentBins`.
- **Slim BCF keep-sets**: SVTYPE-specific — intervals `{END, MATCHA_BOFF}`, BND `{CHR2, POS2, MATCHA_BOFF}`. SVTYPE and SVLEN are not written (SVTYPE is encoded in the filename; SVLEN is derivable and matchcore never reads it). REF, ALT, and QUAL are also blanked (`N`, `.`, missing) since matchcore never reads them. Two slim header templates are built once per run (`buildSlimHdr`) — all FORMAT defs and non-keep INFO defs stripped — and duped per writer.

---

## Preprocessing

Each record is normalized, validated, and slimmed into a temporary `matcha_<pid>_<A|B>_<SVTYPE>_b<bin>.bcf`. Writers are opened lazily. After preproc, each temp BCF is indexed (CSI).

### Field resolution

| Field | Rule | Failure action |
|---|---|---|
| **SVTYPE** | Prefer `INFO/SVTYPE`; fall back to symbolic ALT or bracket BND. **ALT wins on disagreement.** | Unresolvable → warn+skip. `TRA` → warn+skip. `INS` → silent count. |
| **END** (intervals) | Prefer `INFO/END`; fall back to `POS + abs(SVLEN)`. | Neither → warn+skip. `END ≤ POS` → warn+skip. |
| **SVLEN** (intervals) | Prefer `abs(INFO/SVLEN)`; fall back to `END − POS`. | If both provided and disagree >10%, warn and use `END − POS`. |
| **BND mate** | Parse `CHR2`/`POS2` from ALT bracket notation; **ALT is authoritative** — overwrites any stale INFO. `endPos = POS + 1` is used only for `hts_idx_push` (CSI requires endPos > pos); neither END nor SVLEN is written to the slim record. | Malformed ALT → warn+skip. |
| **ID** | Synthesize `CHROM_POS_SVTYPE_LINENUMBER` when absent or `.`. | Always succeeds. |
| **MATCHA_BOFF** | BGZF virtual offset of the source record encoded as `Number=2 Integer` (high32, low32). Points back into the *original* file. | Always succeeds. |
| **Size bin** | `binIndexFor(svlen)`: log2 scale from 1024 bp. BND always lands in bin 0. | Clamped to 0 for non-positive SVLEN. |

Slim BCF content: intervals carry only `END` + `MATCHA_BOFF`; BND carry only `CHR2` + `POS2` + `MATCHA_BOFF`. REF is set to `N`, ALT to `.`, QUAL to BCF missing — all unused by matchcore. FORMAT defs and all other INFO fields are stripped from the writer header.

Warnings go to stderr with prefix `[matcha preproc WARN]`, throttled at 5 per reason per callset (override: `MATCHA_WARN_CAP`).

### Work queue

Jobs are `(chrom, svtype, binA)` tuples. `adjacentBins(binA, threshold, populatedB)` returns the B bins whose size range can overlap under the active threshold — jobs with no adjacent populated B bins are skipped. BND jobs always use `binA = 0` and pair against the single `(svBND, 0)` B BCF.

---

## Size bins

Log2 scale: bin 0 = `[0, 1024)` bp; bin N (N ≥ 1) = `[2^(N+9), 2^(N+10))` bp. See [src/matcha/bins.nim](src/matcha/bins.nim) for `binIndexFor`, `adjacentBins`, and derivations.

---

## Matching

### Interval path (DEL/DUP/INV) — `matchcore.streamJobPairs`

Uses a **TiledBuffer** per B bin — a lazy, eviction-based cache keyed by tile index. Position window per A record: asymmetric `[posA − U, posA + svlenA)` where `U` = upper bound of the B bin's size range. Tiles of width `U` are loaded via CSI region query and evicted once A has advanced past any possible use. Records straddling tile boundaries are assigned to the tile containing their POS (avoids double-counting).

### BND path — `matchcore.streamBndJobPairs`

Maintains a `Deque[BndCacheRec]` of B records in the window `(posA − slop, posA + slop)`. On each A advance: evict left-side records, then CSI-query only the delta `[cacheEnd, posA + slop)` on the right. Each B record is therefore fetched and decoded at most once across all A records whose windows include it.

### Thread pool

Single shared atomic counter; workers `fetchAdd` to claim job indices and write results into disjoint slots. With `--threads 1` the pool is bypassed and jobs run inline. With `--threads ≥ 2`, preproc of A and B runs in two parallel threads before the work queue is built.

### Self-mode dedup (`--self`)

`filesB` aliases `filesA` — single preproc pass. Per-pair filter `aOff < bOff` on `MATCHA_BOFF` eliminates the symmetric `(Y,X)` duplicate and the trivial self-self case. Interval work queue prunes `binsB` to `{b ≥ binA}` so cross-bin pairs are built once; BND path uses the same `aOff < bOff` filter in its emit callback.

---

## Output assembly

`match` mode: main thread iterates per-job result slots in sorted order (CHROM then SVTYPE) and writes TSV rows with `##matcha_metric=` preamble + `#`-header.

`anno` mode: matches are grouped into `Table[aOffset, seq[AnnoMatch]]`. The original input file is reopened and streamed; `bgzf_tell` offset is the join key. Output format auto-detected from `-o` extension (`.vcf`, `.vcf.gz`, `.bcf`; default stdout VCF). Bgzipped outputs get a `.csi` index.

---

## Source-record retrieval (MATCHA_BOFF pattern)

`INFO/MATCHA_BOFF` encodes the BGZF virtual offset as `Number=2 Integer` (high32, low32). Consumers (anno, and future collapse/merge) collect offsets from match results, sort ascending, then `bgzf_seek` + `bcf_read` in order — equivalent amortisation to a region query, strictly better for sparse retrievals.

---

## Module map

| File | Responsibility |
|---|---|
| [src/matcha/main.nim](src/matcha/main.nim) | CLI parsing (`std/parseopt`), subcommand dispatch |
| [src/matcha/utils.nim](src/matcha/utils.nim) | Shared types: `SvType`, `Metric`, `MatchResult`, `MatchConfig` |
| [src/matcha/intervals.nim](src/matcha/intervals.nim) | `reciprocalOverlap`, `jaccard` |
| [src/matcha/bins.nim](src/matcha/bins.nim) | `binIndexFor`, `adjacentBins`, `TiledBuffer` |
| [src/matcha/preproc.nim](src/matcha/preproc.nim) | Normalize → per-(svtype,bin) BCF + work queue; BND ALT parsing; `extraKeepInfo` for anno |
| [src/matcha/matchcore.nim](src/matcha/matchcore.nim) | `streamJobPairs[B,R]` (interval, tiled-buffer) and `streamBndJobPairs[B,R]` (BND, deque+delta) |
| [src/matcha/match.nim](src/matcha/match.nim) | match-mode adapter: build results, self-mode dedup, thread pool, TSV output |
| [src/matcha/anno.nim](src/matcha/anno.nim) | anno-mode: expression parser, `applyAggFunc`, B-side INFO extraction, output VCF assembly |
| [src/matcha/log.nim](src/matcha/log.nim) | Verbose logging (stderr, timestamped) |

---

## Testing

Fixtures in [tests/fixtures/](tests/fixtures/) generated by [tests/generate_fixtures.py](tests/generate_fixtures.py) (needs `bcftools` + `bgzip`). Expected TSVs cover default / strict / jaccard-only / self thresholds.

`nimble test` — each test prints `<ID>\tPASS|FAIL\t<elapsed>\t<desc>`. Per-test timeout via the `timed` template in `test_utils.nim` (default 10s; override with `MATCHA_TEST_TIMEOUT`). Tests mirror module names.
