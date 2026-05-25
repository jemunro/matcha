# matcha — Design & Architecture

Reference for contributors and Claude Code. For user documentation see [README.md](README.md).

---

## Pipeline overview

Three phases per invocation (collapse adds two more):

1. **Preproc** — normalize + slim each input into per-(SVTYPE, bin) temp BCFs, indexed with CSI. Each written record gets a sequential `INFO/SRC_INDEX` (see invariants). `buildWorkQueue` returns `(jobs: seq[MatchJob], fileList: seq[string])`; `fileList` is the deduplicated ordered list of all slim BCF paths used as a `FILE_IDX` lookup table.
2. **Match** — thread pool processes a work queue of `(chrom, svtype, binA)` jobs; each job restricts A and B queries to `job.chrom` using TiledBuffer (intervals), a sliding deque (BND), or a sliding deque (INS). `matchcore.streamJobPairs` / `streamBndJobPairs` / `streamInsJobPairs` return 32-byte `MatchPair` structs — each A record also stamps FILTER, QUAL (Q14.2 uint16), and CALLER_IDX onto the pair, so collapse/merge can build `passQualMap` from MatchPairs alone.
3. **Resolve + Output** (per-mode):
   - `match` → main thread opens slim BCF handles indexed by `fileList`; per-pair CSI `chrom:pos-pos` query retrieves END + ID; writes TSV directly. Header: `#CHROM_A POS_A END_A ID_A CHROM_B POS_B END_B ID_B SVTYPE SIMILARITY`; chrom resolved via `chromOrder[pair.chromIdx]`; `CHROM_A == CHROM_B` always; BND rows emit `.` for END columns.
   - `anno` → B records retrieved per-match via CSI `chrom:posB-posB` query. DB INFO values needed for aggregation are captured during B retrieval from slim B BCFs (which carry user-requested fields via `extraKeepInfo`), so the original DB file is never reopened. Phase 3 streams original A, joining by incrementing SRC_INDEX counter. Output format auto-detected from `-o` extension; bgzipped outputs get a `.csi` index.
   - `collapse` → similarity map from MatchPairs; singletons emitted by matchcore; `passQualMap` built directly from MatchPair `passA`/`qualQ`/`callerIdxA` fields; representative records retrieved via CSI for output.

---

## Key invariants

- **SVTYPE resolution**: ALT wins over `INFO/SVTYPE` on disagreement; ALT is authoritative for BND mate (`CHR2`/`POS2`).
- **SRC_INDEX**: `INFO/SRC_INDEX` (`Number=1 Integer`) is written to every slim BCF record. It is a sequential counter incremented for every record read (including skipped ones), so an equivalent counter loop over the original A file produces matching values in anno's phase-3 join. CSI `chrom:pos-pos` queries with SRC_INDEX as tiebreaker provide O(1) record retrieval from slim BCFs.
- **`--chrs` filter**: Applied before any field resolution (preproc record loop; `integratedMerge` synced-reader loop). Records on excluded chroms are dropped, as are BND records whose mate `CHR2` is excluded. Output `##contig` lines are restricted to kept chroms via `addContigsUnion` (VCF/BCF modes); `anno` additionally strips excluded contigs from the copied input header and skips excluded records in the phase-3 output pass while still advancing the `SRC_INDEX` counter (so the join stays consistent). A post-preproc warning lists any `--chrs` entries absent from all input headers.
- **Self-mode dedup**: `srcIndexA < srcIndexB` filter eliminates symmetric pairs and self-self. The filter is baked into `matchcore` itself (gated on `cfg.selfMode`); no adapter-side work needed. Interval work queue prunes `binsB` to `{b ≥ binA}`.
- **BND and INS always in bin 0**: Both are point events; neither is size-binned into the log2 bins used for intervals. They always land in `(svBND, 0)` / `(svINS, 0)` temp BCFs. `adjacentBins` for bin 0 returns the single bin 0 entry.
- **Slim BCF keep-sets**: SVTYPE-specific — intervals `{END, SRC_INDEX}`, BND `{CHR2, POS2, SRC_INDEX}`, INS `{SVLEN, SRC_INDEX}`. SVTYPE is encoded in the filename. REF, ALT, and QUAL are blanked on match/anno slim BCFs; in collapse/merge paths, BND and INS source REF/ALT are preserved verbatim (`preserveBndAlt`/`preserveInsAlt`) so sequence-resolved insertions and strand-encoded BND bracket forms survive into cohort output. Anno passes `extraKeepInfo` so user-requested DB INFO fields are kept on slim B records, avoiding re-opening the original DB file during resolution.
- **MatchPair contract**: the three matchcore streamers return a flat `seq[MatchPair]` (POD struct, see [src/matcha/utils.nim](src/matcha/utils.nim)). Each pair carries A-side provenance (`passA`, `qualQ`, `callerIdxA`) and `chromIdx` (inherited from `job.chromIdx`, the same chrom for both A and B). This lets collapse/merge build `passQualMap` without a second slim-BCF pass. `qualQ` is Q14.2 fixed-point (`[0, 16383.75]`, precision `0.25`). Resolution is always main-thread.

---

## Preprocessing

Each record is normalized, validated, and slimmed into a temporary `matcha_<pid>_<A|B>_<SVTYPE>_b<bin>.bcf`. Writers are opened lazily. After preproc, each temp BCF is CSI-indexed.

Temp BCFs live inside a per-invocation subdirectory created via `makeRunTmpDir` (`preproc.nim`), which wraps `std/tempfiles.createTempDir` (POSIX `mkdtemp` — atomic, unique even under PID reuse across shared scratch / cluster nodes). All four subcommands create this subdir up front and `removeDir` it on completion. The parent directory is `--tmp-dir` if supplied, otherwise the system temp dir.

### Field resolution

| Field | Rule | Failure action |
|---|---|---|
| **SVTYPE** | Prefer `INFO/SVTYPE`; fall back to symbolic ALT or bracket BND. **ALT wins on disagreement.** | Unresolvable → warn+skip. `TRA` → warn+skip. |
| **END** (intervals) | Prefer `INFO/END`; fall back to `POS + abs(SVLEN)`. | Neither → warn+skip. `END ≤ POS` → warn+skip. |
| **SVLEN** (intervals) | Prefer `abs(INFO/SVLEN)`; fall back to `END − POS`. | If both provided and disagree >10%, warn and use `END − POS`. |
| **BND mate** | Parse `CHR2`/`POS2` from ALT bracket notation; **ALT is authoritative** — overwrites any stale INFO. `endPos = POS + 1` is used only for `hts_idx_push`; neither END nor SVLEN is written to the slim record. | Malformed ALT → warn+skip. |
| **INS length** | Fallback chain: `INFO/INSLEN` → `INFO/SVLEN` → `len(ALT) − len(REF)` (plain-sequence ALTs only, not symbolic/BND) → `len(LEFT_SVINSSEQ) + len(RIGHT_SVINSSEQ)`. Written as `INFO/SVLEN` in the slim record. `endPos = POS + 1` (point event). INS records land in `(svINS, 0)`. | No resolvable length → `skUnresolvableInsLen`, warn+skip. |
| **ID** | Synthesize `CHROM_POS_SVTYPE_LINENUMBER` when absent or `.`. | Always succeeds. |
| **SRC_INDEX** | Sequential int32 counter, incremented per record read (before skip checks). Written to slim BCF INFO as `Number=1 Integer`. | Always succeeds. |
| **Size bin** | `binIndexFor(svlen)`: log2 scale from 1024 bp. BND always lands in bin 0. | Clamped to 0 for non-positive SVLEN. |
| **Chrom filter** | When `--chrs` is supplied, records on non-listed chromosomes are rejected before any other field resolution. BND records whose mate `CHR2` is on an excluded chromosome are also rejected (after BND resolution). Counted as `skChromFiltered`. | Silent skip; reflected in per-callset summary. |

Warnings go to stderr with prefix `[matcha preproc WARN]`, throttled at 5 per reason per callset (override: `MATCHA_WARN_CAP`). After preprocessing completes, a single warning is emitted listing any `--chrs` entries that were not found in any input header (indicating a possible typo).

### Work queue

`buildWorkQueue(a, b, cfg)` returns `(jobs: seq[MatchJob], fileList: seq[string])`. Jobs are `(chrom, svtype, binA)` tuples. Each `MatchJob` carries `chromIdx: int16` (index into chrom order), `fileIdxA: int16`, and `binsB: Table[int, BinBEntry]` where `BinBEntry = (path: string, fileIdx: int16)`. `fileList` is a globally deduplicated ordered list of all slim BCF paths.

Job pruning uses two distinct mechanisms:

- **B-side adjacency**: `adjacentBins(binA, threshold, populatedB)` returns the B bins whose size range can overlap under the active threshold — jobs with no adjacent populated B bins are skipped. BND jobs always use `binA = 0` and pair against the single `(svBND, 0)` B BCF.
- **A-side `(svt, binA, chrom)` triple filter**: `PreprocOutput.populated: Table[SvtypeBin, HashSet[string]]` records exactly which chromosomes have ≥1 record in each `(svtype, bin)` slim BCF. `buildWorkQueue` emits a chrom job only when `chrom ∈ populated[(svt, binA)]`, so chromosomes where the A BCF has no records for that bin are skipped without opening the file. `IntegratedMergeResult.populated` carries the same shape for the collapse/merge paths.

Jobs are sorted **largest A BCF first** (`getFileSize(job.pathA)` descending) as an LPT (longest-processing-time) heuristic for better thread-pool load balance, then **longer chrom first** (`PreprocOutput.chromLens`), then `(svtype, binA)` for determinism.

---

## Size bins

Log2 scale: bin 0 = `[0, 1024)` bp; bin N (N ≥ 1) = `[2^(N+9), 2^(N+10))` bp. See [src/matcha/bins.nim](src/matcha/bins.nim) for `binIndexFor`, `adjacentBins`, and derivations.

---

## Matching

### Interval path (DEL/DUP/INV) — `matchcore.streamJobPairs`

Uses a **TiledBuffer** per B bin — a lazy, eviction-based cache keyed by tile index. Position window per A record: asymmetric `[posA − U, posA + svlenA)` where `U` = upper bound of the B bin's size range. Tiles of width `U` are loaded via CSI region query and evicted once A has advanced past any possible use. Records straddling tile boundaries are assigned to the tile containing their POS (avoids double-counting). Each passing pair is appended to the result `seq[MatchPair]` directly.

### BND path — `matchcore.streamBndJobPairs`

Maintains a `Deque[BndCacheRec]` of B records in the window `(posA − slop, posA + slop)`. On each A advance: evict left-side records, then CSI-query only the delta `[cacheEnd, posA + slop)` on the right. Each B record is fetched and decoded at most once across all A records whose windows include it.

### INS path — `matchcore.streamInsJobPairs`

Identical sliding-cache structure as BND (both use the shared `advanceSlidingCache[T]` generic helper). Cache records carry `pos` + `svlen`. For each A/B candidate pair: `pos_sim = (slop − |dPOS|) / slop`; `len_sim = min(svlenA, svlenB) / max(svlenA, svlenB)`; `sim = sqrt(pos_sim · len_sim)`. Pairs with `sim < cfg.insMinSim` are rejected. Self-mode dedup is the same `srcIndexA < srcIndexB` filter as BND.

### Thread pool

Single shared atomic counter; workers `fetchAdd` to claim job indices and write results into disjoint slots. With `--threads 1` the pool is bypassed and jobs run inline. With `--threads ≥ 2`, preproc of A and B runs in two parallel threads before the work queue is built.

`runMatchPairJobsWithPool` — workers return `seq[seq[MatchPair]]` (pair-only, no resolution). Used by all three modes. Resolution is always handled by the calling mode on the main thread.

### Self-mode dedup (`--self`)

`filesB` aliases `filesA` — single preproc pass. Dedup itself is covered by the `srcIndexA < srcIndexB` invariant above. Collapse always runs in self mode with `emitSingletons = true`.

---

## Collapse pipeline

`matcha collapse` runs over N single-sample caller VCF/BCFs and produces one representative record per cluster. Steps (see [src/matcha/collapse.nim](src/matcha/collapse.nim)):

1. **`resolveHeaders`** (`mergecore.nim`) — analyse N input headers, build a `MergedHeader` driving INFO/FORMAT defs in the output. Per field:
   - **Compatible** (same Number + same Type across every caller that defines the field): single merged def is emitted; no renames.
   - **Number-only conflict** (same Type, different Number — e.g. one caller's `INFO/SVLEN` is `Number=1` and another's is `Number=.`): silently widened to `Number=.`. No rename. Applies symmetrically to INFO and FORMAT.
   - **Type conflict**: each caller's instance is renamed to `FIELD_<callerName>` in the merged BCFs and the output; the slim-BCF rewrite happens during `integratedMerge` via `applyInfoRename`/`applyFmtRename`. Per-caller renamed defs are added to the output header.
   - **Warning gating**: `resolveHeaders` accepts the `--info` / `--format` output filter lists. Conflict warnings (INFO Type, FORMAT Number, FORMAT Type) fire only when the field would appear in the final output filter; the `headerLines` and rename tables are still populated unconditionally so slim BCFs stay correct for fields the user filtered out.
2. **`integratedMerge`** — fused preproc+merge in one `synced_bcf_reader` (`bcf_srs_t`) pass: streams all N callers in lockstep, normalizes each record, applies INFO/FORMAT renames (conflict resolution), filters to user-selected fields, writes per-(svtype, bin) merged slim BCFs. Assigns `SRC_INDEX` (global sequential counter across all callers) and `CALLER_IDX` (0-based caller index) as INFO fields. Uses a shared `htsThreadPool` for parallel BGZF I/O. When `--chrs` is set, records on excluded chromosomes are dropped here (before normalization); BND records whose `CHR2` is on an excluded chromosome are also dropped.
3. **Pass 1 — self-match** — `runMatchPairJobsWithPool` (self-mode, `emitSingletons = true`) over the merged slim BCFs. Singletons (unmatched records) are included in the MatchPair list, so every record is accounted for without a separate enumeration pass.
4. **`buildSimilarityMap` + `locByIdx`** — `buildSimilarityMap` builds `Table[(int32,int32), float64]` (singletons skipped). `locByIdx: Table[int32, (chromIdx, pos, fileIdx)]` is built from all MatchPairs for later slim-BCF retrieval.
5. **Clustering** — union-find over the similarity map yields components; each component is then agglomeratively clustered (`--linkage` average/single/complete) at the threshold.
6. **`passQualMap`** — built by folding over the MatchPair list: each pair carries `passA`, `qualQ`, and `callerIdxA` stamped by matchcore. Because `emitSingletons = true`, every offset appears as `srcIndexA` in at least one pair, so no second slim-BCF pass is needed.
7. **Representative selection** — `selectRepresentative` walks the `--priority` cascade (`PASS, QUAL, CENTRE, ORDER`); `CENTRE` uses pairwise similarities from `simMap` (prefer record with highest mean similarity to other cluster members), and `callerIdx` from `passQualMap` drives `ORDER` (prefer earlier CLI callers). `ORDER` is always appended as the final tiebreaker.
8. **Output** — scan merged slim BCFs, pick representative records by SRC_INDEX, populate `CALLERS` / `N_CALLERS` / `N_MERGED`, drop internal INFO fields, sort by `(chromOrder, POS)`, write VCF/BCF + optional CSI.

Threading: preproc+merge runs in one integrated streaming pass (`integratedMerge`) with a shared `htsThreadPool`; matching reuses the same thread pool as `matcha match`.

---

## Agglomerative clustering

Used by both collapse and merge. `agglomerateComponent(offsets, simMap, linkage, threshold)` handles each connected component; `selfMatchAndCluster` orchestrates the full pipeline.

### Connected components

`buildComponents(simMap)` runs union-find over all (i, j) pairs in the similarity map. Returns `byComp: Table[int32, seq[int32]]` — one entry per component root, value is the sorted list of offsets in that component.

Connected components can grow larger than the direct-match size-bin window suggests. Adjacent size bins can produce edges between bin N and bin N+1; a chain of such edges (bin 7 ↔ 8 ↔ … ↔ 18) transitively links the entire size spectrum into one component via union-find, even though the smallest and largest records in the chain would never directly match. This is expected — agglomeration then correctly splits the chained component at the threshold. The problem only arises when caller artefacts produce thousands of records in a single component (see large-component warning below).

### Dispatcher

`agglomerateComponent(offsets, simMap, linkage, threshold)` dispatches based on component size:

- `len(offsets) ≤ AggDenseThreshold` (256) → `agglomerateDense`
- `len(offsets) > 256` → `agglomerateSparse`

Both procs return `seq[seq[int32]]` with the same semantics. `agglomerateDense` is the reference implementation and is exercised by all existing test fixtures.

### Dense path — `agglomerateDense`

O(N²) space (full N×N similarity matrix), O(N³) time (full rescan per merge). Correct for N up to ~200; unacceptable at N ≈ thousands.

Each iteration scans all active pairs for the maximum similarity; if it meets the threshold the two clusters are merged using the appropriate Lance-Williams update:

- **single**: `max(dAX, dBX)`
- **complete**: `min(dAX, dBX)`
- **average**: `(sA·dAX + sB·dBX) / (sA + sB)`

Missing edges in the simMap are treated as 0.0 — absent pairs never trigger a merge.

### Sparse path — `agglomerateSparse`

O((N + E) log N) time, O(E) space (plus heap accumulation). Used when N > 256.

Data structures:
- `neighbors: seq[Table[int32, float64]]` — sparse adjacency indexed by cluster id, initialised from simMap edges within this component only.
- `heap: HeapQueue[SparseHeapEntry]` — Nim min-heap with negated similarity so the largest edge pops first. Entry tuple: `(negSim, i, j, vi, vj)` — tuple lex-order gives a smallest-`(i, j)` tiebreak, matching the dense path's row-major scan tiebreak for determinism.
- `version: seq[int32]` — one per cluster, incremented on every merge. Stale heap entries are detected and skipped on pop when `version[i] != vi` or `version[j] != vj` (lazy invalidation — the heap is never rebuilt).

Initial heap: one entry per simMap pair with `sim ≥ threshold`. Pop loop: skip if either cluster is inactive or versions don't match; otherwise the popped entry is the current global maximum — merge and continue. Merge (j into i, keep i): apply Lance-Williams update for all `x ∈ neighbors[i] ∪ neighbors[j]` (excluding i and j), write updated weight to `neighbors[i][x]` and `neighbors[x][i]`, remove `neighbors[x][j]`, push new heap entry iff updated weight `≥ threshold`. Mark cluster j inactive; append `members[j]` to `members[i]`; `inc version[i]`. Terminate when the heap is empty or the next valid pop is below threshold.

All three linkages are reducible under Lance-Williams, so the heap-based "pop current global max and merge" strategy produces the same dendrogram as the dense rescan — output is bit-identical. Absent edges in the simMap are treated as 0.0, same as the dense path.

### Large-component warning

Emitted in `selfMatchAndCluster` (mergecore.nim) per component with `len(offsets) ≥ LargeComponentWarn` (500), before agglomeration runs on that component:

```
[WARN T+s] <mode>: large cluster component: <chrom>/<svtype> N=<N> dominant=<caller>:<pct>% — possible caller artifact
```

Chrom and svtype are invariant within a connected component (matching is partitioned by both), so they are read from any one member. Caller attribution is computed lazily — a `Table[int32, Meta]` is built from the `allPairs` A-side fields only on the first large component encountered; typical runs with no large components pay zero overhead.

---

## Merge pipeline

`matcha merge` produces a cohort multi-sample pVCF from N single-sample SV
callsets (typically `matcha collapse` outputs). Steps (see
[src/matcha/merge.nim](src/matcha/merge.nim)):

1. **`validateMergeInputs`** — open each input; enforce exactly 1 sample
   column per file; reject duplicate sample IDs. Returns
   `sampleIdByCaller` in CLI order.
2. **GT auto-add** — if `--format` omits `GT`, prepend it (cohort
   AC/AN/AF requires GT).
3. **`resolveHeaders`** — same as collapse: build a `MergedHeader` from
   the N input headers with conflict resolution.
4. **`buildSlimHdr`** — header for the slim BCF writers. 1 dummy sample
   column `SAMPLE`, plus `FORMAT/SID` (`Number=1,Type=String`) carrying
   the source sample's ID, plus user `--format` fields and standard SV
   INFO defs.
5. **`buildOutputHdr`** — header for the final pVCF. N sample columns
   named from `sampleIdByCaller`, plus cohort `AC`/`AN`/`AF` INFO defs,
   plus conditional `CALLERS` / `N_CALLERS` (only when any input header
   declares `INFO/CALLERS`). No `FORMAT/SID` (slim-internal only) and no
   `N_MERGED` (collapse-only).
6. **`integratedMerge`** with `stampSID=true`, `preserveBndAlt=true`, `preserveInsAlt=true`:
   - Each slim record carries `FORMAT/SID = sampleIdByCaller[ci]`.
   - BND records keep their source ALT verbatim (bracket form encodes strand orientation not derivable from CHR2/POS2).
   - INS records keep their source REF/ALT verbatim when present (preserves sequence-resolved insertions). Non-BND/INS records blank REF/ALT to `N,.`.
7. **Self-match → cluster** — identical pipeline to collapse: pair
   matchcore (selfMode + emitSingletons), `buildSimilarityMap`,
   `selfMatchAndCluster`, `passQualMap` from MatchPair fields (no second
   slim-BCF pass), `selectRepresentative`.
8. **`writeMergeOutput`** — per cluster:
   - Retrieve each cluster member's slim record via CSI; read
     `FORMAT/SID` to determine output sample column.
   - Build a fresh `bcf1_t` in `outputHdr` space:
     - CHROM, POS, ID, QUAL, FILTER from representative.
     - REF=`N`; ALT reconstructed (`<DEL>`/`<DUP>`/`<INV>` from SVTYPE;
       BND uses representative's preserved bracket-form ALT).
     - INFO copied from representative through `keepInfoForMergeOut` filter.
     - FORMAT per-field: read each member's `K`-element data; build
       `N × maxK` buffer; copy member's data into its sample column,
       fill missing samples with appropriate sentinels (`bcf_gt_missing=0`
       for GT, `bcf_int32_missing` for other ints, NaN-tagged missing for
       floats, `.` for strings). With `--missing-to-ref` (`MergeConfig.missingToRef`),
       the GT sentinel for absent samples becomes BCF-encoded REF unphased
       (`2` = `(0+1)<<1 | 0`) instead of `bcf_gt_missing`, so the
       downstream AC/AN accumulator naturally counts them as `0/0`. Only
       absent samples (`mi < 0`) are affected; in-call missing alleles
       from present samples still write `0` and are skipped.
     - Cohort INFO: AC/AN/AF computed from the assembled GT array
       (BCF-encoded; `bcf_gt_missing` skipped, otherwise allele =
       `(v shr 1) - 1`).
     - CALLERS union: when any cluster member's slim record carries
       `INFO/CALLERS` (e.g. from a prior `matcha collapse` run), union the
       caller names preserving representative-first order; `N_CALLERS` =
       distinct count.
   - Same-sample cluster collisions (rare, not blocked at clustering)
     are resolved at this stage: when two cluster members share a
     sample ID, the priority cascade picks one for that column and a
     throttled warning is emitted.
   - Sort buffer by `(chromOrderIdx, pos)`; write to output VCF/BCF +
     optional CSI.

---

## Module map

| File | Responsibility |
|---|---|
| [src/matcha/main.nim](src/matcha/main.nim) | CLI parsing, subcommand dispatch |
| [src/matcha/utils.nim](src/matcha/utils.nim) | Shared types (`SvType`, `Metric`, `MatchPair`, `MatchConfig`, `OutputHeader`) and helpers (`quantizeQual`) |
| [src/matcha/intervals.nim](src/matcha/intervals.nim) | `reciprocalOverlap`, `jaccard` |
| [src/matcha/bins.nim](src/matcha/bins.nim) | Size-bin indexing, `TiledBuffer`, `BufferedRec` |
| [src/matcha/preproc.nim](src/matcha/preproc.nim) | Normalize + slim into per-(svtype, bin) temp BCFs; build the matching work queue |
| [src/matcha/matchcore.nim](src/matcha/matchcore.nim) | Three streamers (interval / BND / INS) returning `seq[MatchPair]`; slim-BCF decode helpers |
| [src/matcha/match.nim](src/matcha/match.nim) | match-mode driver: pair-only pool, chr:pos CSI resolution, TSV output |
| [src/matcha/anno.nim](src/matcha/anno.nim) | anno-mode driver: expression parser, aggregation, phase-3 join, VCF assembly |
| [src/matcha/mergecore.nim](src/matcha/mergecore.nim) | Header merge, clustering, representative selection; `selfMatchAndCluster` shared pipeline |
| [src/matcha/collapse.nim](src/matcha/collapse.nim) | collapse-mode driver: `integratedMerge`, `selfMatchAndCluster`, output assembly |
| [src/matcha/merge.nim](src/matcha/merge.nim) | merge-mode driver: cohort pVCF assembly, per-sample FORMAT routing via `FORMAT/SID`, AC/AN/AF computation |
| [src/matcha/log.nim](src/matcha/log.nim) | Stderr logging, `warnCap` throttle |
| [src/matcha/synced_bcf_reader.nim](src/matcha/synced_bcf_reader.nim) | FFI bindings for htslib `bcf_srs_t`; [csrc/synced_bcf_wrap.c](csrc/synced_bcf_wrap.c) macro wrappers |

---

## Memory management

hts-nim's `VCF.close()` only closes the file handle; the BCF header, CSI/tabix indexes, and `bcf1_t` buffer are left for the GC finalizer. Under ORC the finalizer is not guaranteed to run between jobs, so a per-job `open + close` pattern accumulates C-heap — the BCF header alone is ~500 KB per file with realistic INFO/FORMAT lines.

**Invariant**: any hot/repeating path that opens a `VCF` must call `tearDownVcf(v)` (defined in [src/matcha/matchcore.nim](src/matcha/matchcore.nim)) instead of `v.close()`. `tearDownVcf` closes the file and explicitly calls `bcf_hdr_destroy` on the header, dropping ~50× of the per-job leak. Remaining leaks (`bidx`, `tidx`, the `bcf1_t` buffer — private fields of hts-nim's `VCF`) are bounded at ~0.2 MB per opened VCF until the GC eventually runs and are accepted rather than reached into via private layout.

One-shot callsites (preproc, header inspection, output writer) keep using `close()`; the lifetime is the process and the finalizer eventually clears them.

---

## Testing

Fixtures in [tests/fixtures/](tests/fixtures/) generated by [tests/generate_fixtures.py](tests/generate_fixtures.py) (needs `bcftools` + `bgzip`). Expected TSVs cover default / strict / jaccard-only / self thresholds.

`nimble test` — each test prints `<ID>\tPASS|FAIL\t<elapsed>\t<desc>`. Per-test timeout via the `timed` template in `test_utils.nim` (default 10s; override with `MATCHA_TEST_TIMEOUT`). Tests mirror module names.
