# matcha

Compiled SV matching and annotation tool (Nim + hts-nim). See [README.md](README.md) for user docs, [DESIGN.md](DESIGN.md) for architecture detail.

## Status

Milestones 1–4 (`matcha match`, `matcha anno`, `matcha collapse`, `matcha merge`) complete. SVTYPE=INS support is out of scope (silent skip).

## Build & test

```
nimble build          # → ./matcha
nimble test           # run test suite
python3 tests/generate_fixtures.py   # regenerate fixtures (needs bcftools + bgzip)
```

## Module map

| File | Responsibility |
|---|---|
| `main.nim` | CLI (`std/parseopt`), subcommand dispatch |
| `utils.nim` | Shared types: `SvType`, `Metric`, `MatchPair`, `MatchConfig`, `OutputHeader`, `NO_MATCH` |
| `intervals.nim` | `reciprocalOverlap`, `jaccard` |
| `bins.nim` | `binIndexFor`, `adjacentBins`, `TiledBuffer`, `BufferedRec` |
| `preproc.nim` | Normalize → per-(svtype,bin) temp BCF + work queue; `extraKeepInfo` (anno); `SRC_INDEX` assignment; `buildWorkQueue` returns `(jobs, fileList)`; `parseChrsArg`/`warnMissingChrs` for `--chrs` filter |
| `matchcore.nim` | `streamJobPairs` (interval) and `streamBndJobPairs` (BND) — both return `seq[MatchPair]` (28-byte: srcIndex, pos, sim, fileIdx, chromIdx, svtype). Slim-BCF decode helpers (`readSrcIndex`, `readPos2`, `readChr2`, `extractEnd`). |
| `match.nim` | match-mode: pair-only pool (`runMatchPairJobsWithPool`), main-thread chr:pos CSI resolution, TSV output |
| `anno.nim` | anno-mode: expression parser, `applyAggFunc`, per-match chr:pos CSI B retrieval, SRC_INDEX counter join for phase 3, output VCF assembly |
| `mergecore.nim` | header merge (`resolveHeaders`), `buildSimilarityMap`, union-find, agglomerative clustering, `selectRepresentative` |
| `collapse.nim` | collapse-mode driver: `integratedMerge` (fused preproc+merge via `synced_bcf_reader`), self-match with singleton emission, clustering, representative output |
| `merge.nim` | merge-mode driver: cohort pVCF across N single-sample inputs. Builds slimHdr (1 dummy `SAMPLE` + `FORMAT/SID`) and outputHdr (N samples + AC/AN/AF); per-sample FORMAT routing via SID. |
| `log.nim` | Verbose/warn/error logging (stderr, timestamped); `warnCap` throttle |
| `synced_bcf_reader.nim` | FFI bindings for htslib `bcf_srs_t`; `newVariantView`/`setRecView`; `csrc/synced_bcf_wrap.c` macro wrappers |
