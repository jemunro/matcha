# matcha — Design & Architecture

Reference for contributors and Claude Code. For user documentation see [README.md](README.md).

---

## Pipeline overview

Three phases per invocation (collapse adds two more):

1. **Preproc** — normalize + slim each input into per-(SVTYPE, bin) temp BCFs, indexed with CSI.
2. **Match** — thread pool processes a work queue of `(chrom, svtype, binA)` jobs; each job pairs A records against adjacent B bins using TiledBuffer (intervals) or a sliding deque (BND). `matchcore.streamJobPairs` / `streamBndJobPairs` return a minimal `seq[MatchPair]` of `(aOff, bOff, sim)` triples — no field resolution happens here.
3. **Resolve** (per-mode) — each adapter scans the same slim BCFs once more, picking records by `MATCHA_BOFF`, to materialise the fields it needs:
   - `match` → `MatchResult` (`chromA/posA/endA/idA`, `chromB/posB/endB/idB`, `svtype`, `sim`) for TSV output.
   - `anno`  → `AnnoMatch` carrying user-requested DB INFO values; aggregation runs over these grouped by `aOffset`.
   - `collapse` → skips resolution entirely; `buildSimilarityMap` consumes the raw `MatchPair` triples.
4. **Output** — main thread collects per-job results in deterministic (job-sorted) order and writes TSV (`match`), annotated VCF/BCF (`anno`), or representative-record VCF/BCF (`collapse`).

---

## Key invariants

- **SVTYPE resolution**: ALT wins over `INFO/SVTYPE` on disagreement; ALT is authoritative for BND mate (`CHR2`/`POS2`).
- **MATCHA_BOFF**: BGZF virtual offset of the source record, stored as `Number=2 Integer` (high32, low32). Doubles as a stable record identity used (a) as the join key during per-mode resolution against the slim BCFs, (b) by anno's phase 3 to `bgzf_seek` back into the original input, and (c) by collapse for canonical pair keys and clustering.
- **Self-mode dedup**: `aOff < bOff` filter on `MATCHA_BOFF` eliminates symmetric pairs and self-self. The filter is baked into `matchcore` itself (gated on `cfg.selfMode`); no adapter-side work needed. Interval work queue prunes `binsB` to `{b ≥ binA}`.
- **BND always in bin 0**: BND records are point events, not size-binned. They always land in `(svBND, 0)` temp BCFs and bypass `adjacentBins`.
- **Slim BCF keep-sets**: SVTYPE-specific — intervals `{END, MATCHA_BOFF}`, BND `{CHR2, POS2, MATCHA_BOFF}`. Anno passes `extraKeepInfo` to preprocess the database with the user-requested DB INFO fields preserved on the slim B records (so DB resolution never has to re-open the original). SVTYPE and SVLEN are not written (SVTYPE is encoded in the filename; SVLEN is derivable and matchcore never reads it). REF, ALT, and QUAL are also blanked (`N`, `.`, missing) since matchcore never reads them. Two slim header templates are built once per run (`buildSlimHdr`) — all FORMAT defs and non-keep INFO defs stripped — and duped per writer.
- **Minimal matchcore contract**: `streamJobPairs` and `streamBndJobPairs` return `seq[MatchPair]` (24-byte POD: `aOff`, `bOff`, `sim`). No generic parameters, no `extract`/`emit` callbacks, no per-mode payload cache. Per-mode logic lives in `match.nim` and `anno.nim` resolvers.

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

Uses a **TiledBuffer** per B bin — a lazy, eviction-based cache keyed by tile index. Position window per A record: asymmetric `[posA − U, posA + svlenA)` where `U` = upper bound of the B bin's size range. Tiles of width `U` are loaded via CSI region query and evicted once A has advanced past any possible use. Records straddling tile boundaries are assigned to the tile containing their POS (avoids double-counting). Each passing pair is appended to the result `seq[MatchPair]` directly — no payload-cache table, no callback indirection.

### BND path — `matchcore.streamBndJobPairs`

Maintains a `Deque[BndCacheRec]` of B records in the window `(posA − slop, posA + slop)`. On each A advance: evict left-side records, then CSI-query only the delta `[cacheEnd, posA + slop)` on the right. Each B record is therefore fetched and decoded at most once across all A records whose windows include it.

### Thread pool

Single shared atomic counter; workers `fetchAdd` to claim job indices and write results into disjoint slots. With `--threads 1` the pool is bypassed and jobs run inline. With `--threads ≥ 2`, preproc of A and B runs in two parallel threads before the work queue is built.

Two pool entry points share the same machinery:

- `runMatchJobsWithPool` — workers run match + per-job resolution, returning `seq[seq[MatchResult]]`. Used by `matcha match`.
- `runMatchPairJobsWithPool` — workers return `seq[seq[MatchPair]]` and skip resolution entirely. Used by `matcha collapse` (and any future mode that only needs the pair triples).

### Self-mode dedup (`--self`)

`filesB` aliases `filesA` — single preproc pass. The `aOff < bOff` filter on `MATCHA_BOFF` lives inside `matchcore` itself (interval and BND paths both honour `cfg.selfMode`), eliminating the symmetric `(Y,X)` duplicate and the trivial self-self case. Interval work queue additionally prunes `binsB` to `{b ≥ binA}` so cross-bin pairs are built once. Collapse always runs in self mode.

---

## Output assembly

`match` mode: main thread iterates per-job result slots in deterministic (job-sorted) order (header chrom order, then SVTYPE, then `binA`) and writes TSV rows with `##matcha_metric=` preamble + `#`-header. The header is `#CHROM_A POS_A END_A ID_A CHROM_B POS_B END_B ID_B SVTYPE SIMILARITY`; `CHROM_A == CHROM_B` always under the current per-chrom job model.

`anno` mode: matches are grouped into `Table[aOffset, seq[AnnoMatch]]`. The original input file is reopened and streamed; `bgzf_tell` offset on each read is the join key. DB INFO values needed for aggregation were captured during per-job resolution from the slim B BCFs (which already carry the user-requested fields via `extraKeepInfo`), so the original DB file is never reopened. Output format auto-detected from `-o` extension (`.vcf`, `.vcf.gz`, `.bcf`; default stdout VCF). Bgzipped outputs get a `.csi` index.

`collapse` mode: pair-only match → similarity map → union-find components → per-component agglomerative clustering → representative selection (`--priority` cascade) → output assembly streams the *original* caller files, reads each representative record by `bgzf_seek` + `bcf_read`, sorts by `(chromOrder, POS)`, and writes the merged VCF/BCF with `SOURCE`/`SOURCELIST`/`N_SOURCE`/`N_MERGED` provenance INFO populated.

---

## Source-record retrieval (MATCHA_BOFF pattern)

`INFO/MATCHA_BOFF` encodes the BGZF virtual offset as `Number=2 Integer` (high32, low32). It plays two distinct roles:

- **Slim-BCF join key** — during per-mode resolution, adapters re-scan the per-job slim BCFs and pick records by `MATCHA_BOFF`. Slim BCFs are tiny (keep-set only) so a sequential scan is cheap and avoids any random I/O on the originals.
- **Original-file seek key** — anno's phase 3 streams the original input and uses `bgzf_tell` per record as the join against `Table[aOffset, ...]`; collapse seeks the original caller files (`bgzf_seek` + `bcf_read`) to materialise representative records for output. Collected offsets are sorted ascending before seeking — amortised cost equivalent to a region query, strictly better for sparse retrievals.

---

## Collapse pipeline

`matcha collapse` runs over N single-sample caller VCF/BCFs and produces one representative record per cluster. Steps (see [src/matcha/collapse.nim](src/matcha/collapse.nim)):

1. **`resolveHeaders`** (`mergecore.nim`) — analyse N input headers, build a `MergedHeader` with conflict resolution for INFO/FORMAT defs that disagree across callers.
2. **`integratedMerge`** — fused preproc+merge in one `synced_bcf_reader` (`bcf_srs_t`) pass: streams all N callers in lockstep, normalizes each record, applies INFO/FORMAT renames (conflict resolution), filters to user-selected fields, writes per-(svtype, bin) merged slim BCFs. Uses a shared `htsThreadPool` for parallel BGZF I/O across all readers and writers. MATCHA_BOFF stores a composite (callerIdx, bgzfOffset) token.
3. **Pass 1 — self-match** — `runMatchPairJobsWithPool` (self-mode, no resolution) over the merged slim BCFs returns `seq[MatchPair]` directly. `buildSimilarityMap` builds a canonical `Table[(aOff, bOff), sim]`.
4. **Pass 2 — `exploreMerged`** — sequential scan of the merged slim BCFs to enumerate every offset and capture PASS/QUAL for representative selection.
5. **Clustering** — union-find over the similarity map yields components; each component is then agglomeratively clustered (`--linkage` average/single/complete) at the threshold.
6. **Representative selection** — `selectRepresentative` walks the `--priority` cascade (`PASS, QUAL, CENTRE, ORDER`); `ORDER` is always appended as the final tiebreaker.
7. **Output** — re-open each original caller file, `bgzf_seek` to representative offsets, translate headers to the output, populate `SOURCE` / `SOURCELIST` / `N_SOURCE` / `N_MERGED`, sort by `(chromOrder, POS)`, write VCF/BCF + optional CSI.

Threading: preproc+merge runs in one integrated streaming pass (`integratedMerge`) with a shared `htsThreadPool`; matching reuses the same thread pool as `matcha match`.

---

## Module map

| File | Responsibility |
|---|---|
| [src/matcha/main.nim](src/matcha/main.nim) | CLI parsing (`std/parseopt`), subcommand dispatch |
| [src/matcha/utils.nim](src/matcha/utils.nim) | Shared types: `SvType`, `Metric`, `MatchResult`, `MatchConfig`, `OutputHeader`, `formatMatchResult` |
| [src/matcha/intervals.nim](src/matcha/intervals.nim) | `reciprocalOverlap`, `jaccard` |
| [src/matcha/bins.nim](src/matcha/bins.nim) | `binIndexFor`, `adjacentBins`, `TiledBuffer`, `BufferedRec` |
| [src/matcha/preproc.nim](src/matcha/preproc.nim) | Normalize → per-(svtype,bin) BCF + work queue; BND ALT parsing; `extraKeepInfo` (anno); `MATCHA_BOFF` encoding |
| [src/matcha/matchcore.nim](src/matcha/matchcore.nim) | `streamJobPairs` (interval, tiled-buffer) and `streamBndJobPairs` (BND, deque+delta), both returning `seq[MatchPair]`. Slim-BCF INFO decode helpers (`readBoff`, `readPos2`, `readChr2`, `extractEnd`). |
| [src/matcha/match.nim](src/matcha/match.nim) | match-mode adapter: slim-BCF resolution → `MatchResult`, thread pools (`runMatchJobsWithPool`, `runMatchPairJobsWithPool`), TSV output |
| [src/matcha/anno.nim](src/matcha/anno.nim) | anno-mode: expression parser, `applyAggFunc`, slim-B DB INFO extraction, output VCF assembly |
| [src/matcha/mergecore.nim](src/matcha/mergecore.nim) | Header merge (`resolveHeaders`), `buildSimilarityMap`, union-find components, agglomerative clustering, `selectRepresentative` |
| [src/matcha/collapse.nim](src/matcha/collapse.nim) | collapse-mode driver: `integratedMerge` (fused preproc+merge via `synced_bcf_reader`), Pass 1 matching, Pass 2 exploration, clustering, output assembly |
| [src/matcha/log.nim](src/matcha/log.nim) | Verbose/warn/error logging (stderr, timestamped); `warnCap` throttle |
| [src/matcha/synced_bcf_reader.nim](src/matcha/synced_bcf_reader.nim) | FFI bindings for htslib `bcf_srs_t`; `newVariantView`/`setRecView`; [csrc/synced_bcf_wrap.c](csrc/synced_bcf_wrap.c) macro wrappers |

---

## Testing

Fixtures in [tests/fixtures/](tests/fixtures/) generated by [tests/generate_fixtures.py](tests/generate_fixtures.py) (needs `bcftools` + `bgzip`). Expected TSVs cover default / strict / jaccard-only / self thresholds.

`nimble test` — each test prints `<ID>\tPASS|FAIL\t<elapsed>\t<desc>`. Per-test timeout via the `timed` template in `test_utils.nim` (default 10s; override with `MATCHA_TEST_TIMEOUT`). Tests mirror module names.
