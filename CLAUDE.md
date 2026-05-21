# matcha

Compiled SV matching and annotation tool (Nim + hts-nim). See [README.md](README.md) for user docs, [DESIGN.md](DESIGN.md) for architecture detail.

## Status

Milestones 1–4 (`matcha match`, `matcha anno`, `matcha collapse`, `matcha merge`) complete. SVTYPE=INS is supported across all four subcommands (position + size similarity; `--min-ins-sim`, `--ins-slop`). SVTYPE=TRA is not supported (warned and skipped).

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
| `preproc.nim` | Normalize → per-(svtype,bin) temp BCF + work queue; INS length resolution (`readInsLen`; INSLEN → SVLEN → ALT seq → SVINSSEQ chain); `extraKeepInfo` (anno); `SRC_INDEX` assignment; `buildWorkQueue` returns `(jobs, fileList)`; `parseChrsArg`/`warnMissingChrs` for `--chrs` filter |
| `matchcore.nim` | `streamJobPairs` (interval), `streamBndJobPairs` (BND), `streamInsJobPairs` (INS) — all return `seq[MatchPair]`. Shared `advanceSlidingCache[T]` helper used by BND and INS. Slim-BCF decode helpers (`readSrcIndex`, `readPos2`, `readChr2`, `extractEnd`, `readSvlen`). |
| `match.nim` | match-mode: pair-only pool (`runMatchPairJobsWithPool`), main-thread chr:pos CSI resolution, TSV output |
| `anno.nim` | anno-mode: expression parser, `applyAggFunc`, per-match chr:pos CSI B retrieval, SRC_INDEX counter join for phase 3, output VCF assembly |
| `mergecore.nim` | header merge (`resolveHeaders`), `buildSimilarityMap`, union-find, agglomerative clustering, `selectRepresentative` |
| `collapse.nim` | collapse-mode driver: `integratedMerge` (fused preproc+merge via `synced_bcf_reader`), self-match with singleton emission, clustering, representative output |
| `merge.nim` | merge-mode driver: cohort pVCF across N single-sample inputs. Builds slimHdr (1 dummy `SAMPLE` + `FORMAT/SID`) and outputHdr (N samples + AC/AN/AF); per-sample FORMAT routing via SID. `--missing-to-ref` writes absent samples as `0/0` (encoded GT = `2`) so AC/AN counts them. |
| `log.nim` | Verbose/warn/error logging (stderr, timestamped); `warnCap` throttle |
| `synced_bcf_reader.nim` | FFI bindings for htslib `bcf_srs_t`; `newVariantView`/`setRecView`; `csrc/synced_bcf_wrap.c` macro wrappers |
