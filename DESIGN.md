# matcha — Design & Architecture

Reference for contributors and Claude Code. For user documentation see [README.md](README.md).

---

## Pipeline overview

Three phases per invocation (collapse adds two more):

1. **Preproc** — normalize + slim each input into per-(SVTYPE, bin) temp BCFs, indexed with CSI. Each written record gets `INFO/SRC_INDEX` (sequential int32, incremented for every record read including skipped ones). `buildWorkQueue` returns `(jobs: seq[MatchJob], fileList: seq[string])`; `fileList` is the deduplicated ordered list of all slim BCF paths used as a `FILE_IDX` lookup table.
2. **Match** — thread pool processes a work queue of `(chrom, svtype, binA)` jobs; each job pairs A records against adjacent B bins using TiledBuffer (intervals) or a sliding deque (BND). `matchcore.streamJobPairs` / `streamBndJobPairs` return 28-byte `MatchPair` triples — no field resolution happens here.
3. **Resolve + Output** (per-mode):
   - `match` → main thread opens slim BCF handles indexed by `fileList`; per-pair CSI `chrom:pos-pos` query retrieves END + ID; writes TSV directly.
   - `anno` → B records retrieved per-match via CSI `chrom:posB-posB` query. Phase 3 streams original A, joining by incrementing SRC_INDEX counter.
   - `collapse` → similarity map from MatchPairs; singletons emitted by matchcore; PASS/QUAL/CALLER_IDX retrieved via targeted CSI queries for cluster members; representative records retrieved the same way for output.

---

## Key invariants

- **SVTYPE resolution**: ALT wins over `INFO/SVTYPE` on disagreement; ALT is authoritative for BND mate (`CHR2`/`POS2`).
- **SRC_INDEX**: `INFO/SRC_INDEX` (`Number=1 Integer`) is written to every slim BCF record. It is a sequential counter incremented for every record read (including skipped ones), so an equivalent counter loop over the original A file produces matching values in anno's phase-3 join. CSI `chrom:pos-pos` queries with SRC_INDEX as tiebreaker provide O(1) record retrieval from slim BCFs.
- **Self-mode dedup**: `srcIndexA < srcIndexB` filter eliminates symmetric pairs and self-self. The filter is baked into `matchcore` itself (gated on `cfg.selfMode`); no adapter-side work needed. Interval work queue prunes `binsB` to `{b ≥ binA}`.
- **BND always in bin 0**: BND records are point events, not size-binned. They always land in `(svBND, 0)` temp BCFs and bypass `adjacentBins`.
- **Slim BCF keep-sets**: SVTYPE-specific — intervals `{END, SRC_INDEX}`, BND `{CHR2, POS2, SRC_INDEX}`. Anno passes `extraKeepInfo` to preprocess the database with user-requested DB INFO fields preserved on slim B records (so DB resolution never has to re-open the original). SVTYPE and SVLEN are not written (SVTYPE is encoded in the filename; SVLEN is derivable). REF, ALT, and QUAL are blanked (`N`, `.`, missing). Two slim header templates are built once per run (`buildSlimHdr`) — all FORMAT defs and non-keep INFO defs stripped — and duped per writer.
- **MatchPair contract**: `streamJobPairs` and `streamBndJobPairs` return 28-byte POD structs `(srcIndexA, srcIndexB: int32; posA, posB: int32; sim: float32; fileIdxA, fileIdxB, chromIdx: int16; svtype: int8)`. No generic parameters, no callbacks, no per-mode payload cache. Resolution is always main-thread.

---

## Preprocessing

Each record is normalized, validated, and slimmed into a temporary `matcha_<pid>_<A|B>_<SVTYPE>_b<bin>.bcf`. Writers are opened lazily. After preproc, each temp BCF is CSI-indexed.

### Field resolution

| Field | Rule | Failure action |
|---|---|---|
| **SVTYPE** | Prefer `INFO/SVTYPE`; fall back to symbolic ALT or bracket BND. **ALT wins on disagreement.** | Unresolvable → warn+skip. `TRA` → warn+skip. `INS` → silent count. |
| **END** (intervals) | Prefer `INFO/END`; fall back to `POS + abs(SVLEN)`. | Neither → warn+skip. `END ≤ POS` → warn+skip. |
| **SVLEN** (intervals) | Prefer `abs(INFO/SVLEN)`; fall back to `END − POS`. | If both provided and disagree >10%, warn and use `END − POS`. |
| **BND mate** | Parse `CHR2`/`POS2` from ALT bracket notation; **ALT is authoritative** — overwrites any stale INFO. `endPos = POS + 1` is used only for `hts_idx_push`; neither END nor SVLEN is written to the slim record. | Malformed ALT → warn+skip. |
| **ID** | Synthesize `CHROM_POS_SVTYPE_LINENUMBER` when absent or `.`. | Always succeeds. |
| **SRC_INDEX** | Sequential int32 counter, incremented per record read (before skip checks). Written to slim BCF INFO as `Number=1 Integer`. | Always succeeds. |
| **Size bin** | `binIndexFor(svlen)`: log2 scale from 1024 bp. BND always lands in bin 0. | Clamped to 0 for non-positive SVLEN. |

Slim BCF content: intervals carry only `END` + `SRC_INDEX`; BND carry only `CHR2` + `POS2` + `SRC_INDEX`. REF is set to `N`, ALT to `.`, QUAL to BCF missing — all unused by matchcore.

Warnings go to stderr with prefix `[matcha preproc WARN]`, throttled at 5 per reason per callset (override: `MATCHA_WARN_CAP`).

### Work queue

`buildWorkQueue(a, b, cfg)` returns `(jobs: seq[MatchJob], fileList: seq[string])`. Jobs are `(chrom, svtype, binA)` tuples. Each `MatchJob` carries `chromIdx: int16` (index into chrom order), `fileIdxA: int16`, and `binsB: Table[int, BinBEntry]` where `BinBEntry = (path: string, fileIdx: int16)`. `fileList` is a globally deduplicated ordered list of all slim BCF paths.

`adjacentBins(binA, threshold, populatedB)` returns the B bins whose size range can overlap under the active threshold — jobs with no adjacent populated B bins are skipped. BND jobs always use `binA = 0` and pair against the single `(svBND, 0)` B BCF.

---

## Size bins

Log2 scale: bin 0 = `[0, 1024)` bp; bin N (N ≥ 1) = `[2^(N+9), 2^(N+10))` bp. See [src/matcha/bins.nim](src/matcha/bins.nim) for `binIndexFor`, `adjacentBins`, and derivations.

---

## Matching

### Interval path (DEL/DUP/INV) — `matchcore.streamJobPairs`

Uses a **TiledBuffer** per B bin — a lazy, eviction-based cache keyed by tile index. Position window per A record: asymmetric `[posA − U, posA + svlenA)` where `U` = upper bound of the B bin's size range. Tiles of width `U` are loaded via CSI region query and evicted once A has advanced past any possible use. Records straddling tile boundaries are assigned to the tile containing their POS (avoids double-counting). Each passing pair is appended to the result `seq[MatchPair]` directly.

### BND path — `matchcore.streamBndJobPairs`

Maintains a `Deque[BndCacheRec]` of B records in the window `(posA − slop, posA + slop)`. On each A advance: evict left-side records, then CSI-query only the delta `[cacheEnd, posA + slop)` on the right. Each B record is fetched and decoded at most once across all A records whose windows include it.

### Thread pool

Single shared atomic counter; workers `fetchAdd` to claim job indices and write results into disjoint slots. With `--threads 1` the pool is bypassed and jobs run inline. With `--threads ≥ 2`, preproc of A and B runs in two parallel threads before the work queue is built.

`runMatchPairJobsWithPool` — workers return `seq[seq[MatchPair]]` (pair-only, no resolution). Used by all three modes. Resolution is always handled by the calling mode on the main thread.

### Self-mode dedup (`--self`)

`filesB` aliases `filesA` — single preproc pass. The `srcIndexA < srcIndexB` filter lives inside `matchcore` itself (interval and BND paths both honour `cfg.selfMode`), eliminating the symmetric `(Y,X)` duplicate and the trivial self-self case. Interval work queue additionally prunes `binsB` to `{b ≥ binA}` so cross-bin pairs are built once. Collapse always runs in self mode with `emitSingletons = true`.

---

## Output assembly

`match` mode: main thread opens one slim BCF handle per file in `fileList`. For each pair, a CSI `chrom:pos-pos` query on the appropriate slim BCF retrieves END and ID; chrom comes from `chromOrder[pair.chromIdx]`, svtype from `SvType(pair.svtype)`. TSV rows are written directly. Header: `#CHROM_A POS_A END_A ID_A CHROM_B POS_B END_B ID_B SVTYPE SIMILARITY`; `CHROM_A == CHROM_B` always; BND rows emit `.` for END columns.

`anno` mode: per-job B retrieval uses targeted CSI `chrom:posB-posB` queries (one seek per unique matched B position, efficient for sparse A vs dense B). Phase 3 streams original A; an incrementing `srcIdx: int32` counter is the join key against `Table[int32, seq[AnnoMatch]]`. DB INFO values needed for aggregation are captured during B retrieval from slim B BCFs (which carry user-requested fields via `extraKeepInfo`), so the original DB file is never reopened. Output format auto-detected from `-o` extension. Bgzipped outputs get a `.csi` index.

`collapse` mode: pair-only self-match (with singleton emission) → similarity map → `locByIdx` built from MatchPairs → clustering → targeted CSI queries for cluster members → representative selection → output scans merged slim BCFs, identifies representatives by SRC_INDEX, strips internal INFO fields, sorts by `(chromOrder, POS)`, writes VCF/BCF.

---

## Collapse pipeline

`matcha collapse` runs over N single-sample caller VCF/BCFs and produces one representative record per cluster. Steps (see [src/matcha/collapse.nim](src/matcha/collapse.nim)):

1. **`resolveHeaders`** (`mergecore.nim`) — analyse N input headers, build a `MergedHeader` with conflict resolution for INFO/FORMAT defs that disagree across callers.
2. **`integratedMerge`** — fused preproc+merge in one `synced_bcf_reader` (`bcf_srs_t`) pass: streams all N callers in lockstep, normalizes each record, applies INFO/FORMAT renames (conflict resolution), filters to user-selected fields, writes per-(svtype, bin) merged slim BCFs. Assigns `SRC_INDEX` (global sequential counter across all callers) and `CALLER_IDX` (0-based caller index) as INFO fields. Uses a shared `htsThreadPool` for parallel BGZF I/O.
3. **Pass 1 — self-match** — `runMatchPairJobsWithPool` (self-mode, `emitSingletons = true`) over the merged slim BCFs. Singletons (unmatched records) are included in the MatchPair list, so every record is accounted for without a separate enumeration pass.
4. **`buildSimilarityMap` + `locByIdx`** — `buildSimilarityMap` builds `Table[(int32,int32), float64]` (singletons skipped). `locByIdx: Table[int32, (chromIdx, pos, fileIdx)]` is built from all MatchPairs for later slim-BCF retrieval.
5. **Clustering** — union-find over the similarity map yields components; each component is then agglomeratively clustered (`--linkage` average/single/complete) at the threshold.
6. **`passQualMap`** — targeted CSI queries (grouped by `fileIdx`) for all cluster members retrieve PASS/QUAL/CALLER_IDX from merged slim BCFs into `Table[int32, (hasPASS, qual, callerIdx)]`.
7. **Representative selection** — `selectRepresentative` walks the `--priority` cascade (`PASS, QUAL, CENTRE, ORDER`); `callerIdx` from `passQualMap` drives `CENTRE` (prefer earlier CLI callers) and `ORDER`. `ORDER` is always appended as the final tiebreaker.
8. **Output** — scan merged slim BCFs, pick representative records by SRC_INDEX, populate `SOURCE` / `SOURCELIST` / `N_SOURCE` / `N_MERGED`, drop internal INFO fields, sort by `(chromOrder, POS)`, write VCF/BCF + optional CSI.

Threading: preproc+merge runs in one integrated streaming pass (`integratedMerge`) with a shared `htsThreadPool`; matching reuses the same thread pool as `matcha match`.

---

## Module map

| File | Responsibility |
|---|---|
| [src/matcha/main.nim](src/matcha/main.nim) | CLI parsing (`std/parseopt`), subcommand dispatch |
| [src/matcha/utils.nim](src/matcha/utils.nim) | Shared types: `SvType`, `Metric`, `MatchPair` (28-byte POD), `MatchConfig`, `OutputHeader`, `NO_MATCH` |
| [src/matcha/intervals.nim](src/matcha/intervals.nim) | `reciprocalOverlap`, `jaccard` |
| [src/matcha/bins.nim](src/matcha/bins.nim) | `binIndexFor`, `adjacentBins`, `TiledBuffer`, `BufferedRec` |
| [src/matcha/preproc.nim](src/matcha/preproc.nim) | Normalize → per-(svtype,bin) BCF + work queue; BND ALT parsing; `extraKeepInfo` (anno); `SRC_INDEX` assignment; `buildWorkQueue` → `(jobs, fileList)` |
| [src/matcha/matchcore.nim](src/matcha/matchcore.nim) | `streamJobPairs` (interval, tiled-buffer) and `streamBndJobPairs` (BND, deque+delta), both returning `seq[MatchPair]`. Slim-BCF INFO decode helpers (`readSrcIndex`, `readPos2`, `readChr2`, `extractEnd`). |
| [src/matcha/match.nim](src/matcha/match.nim) | match-mode: `runMatchPairJobsWithPool`, main-thread chr:pos CSI resolution via `fileList`, TSV output |
| [src/matcha/anno.nim](src/matcha/anno.nim) | anno-mode: expression parser, `applyAggFunc`, per-match chr:pos CSI B retrieval, SRC_INDEX counter phase-3 join, output VCF assembly |
| [src/matcha/mergecore.nim](src/matcha/mergecore.nim) | Header merge (`resolveHeaders`), `buildSimilarityMap`, union-find components, agglomerative clustering, `selectRepresentative` |
| [src/matcha/collapse.nim](src/matcha/collapse.nim) | collapse-mode driver: `integratedMerge` (fused preproc+merge via `synced_bcf_reader`), self-match with singleton emission, `locByIdx`/`passQualMap` via CSI, clustering, output assembly |
| [src/matcha/log.nim](src/matcha/log.nim) | Verbose/warn/error logging (stderr, timestamped); `warnCap` throttle |
| [src/matcha/synced_bcf_reader.nim](src/matcha/synced_bcf_reader.nim) | FFI bindings for htslib `bcf_srs_t`; `newVariantView`/`setRecView`; [csrc/synced_bcf_wrap.c](csrc/synced_bcf_wrap.c) macro wrappers |

---

## Testing

Fixtures in [tests/fixtures/](tests/fixtures/) generated by [tests/generate_fixtures.py](tests/generate_fixtures.py) (needs `bcftools` + `bgzip`). Expected TSVs cover default / strict / jaccard-only / self thresholds.

`nimble test` — each test prints `<ID>\tPASS|FAIL\t<elapsed>\t<desc>`. Per-test timeout via the `timed` template in `test_utils.nim` (default 10s; override with `MATCHA_TEST_TIMEOUT`). Tests mirror module names.
