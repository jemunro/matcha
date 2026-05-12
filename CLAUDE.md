# matcha

Compiled SV matching and annotation tool (Nim + hts-nim). See [README.md](README.md) for user docs, [DESIGN.md](DESIGN.md) for architecture detail.

## Status

Milestones 1–3 (`matcha match`, `matcha anno`, `matcha collapse`) complete. SVTYPE=INS support is out of scope (silent skip). `matcha merge` (cross-sample cohort pVCF) is planned. Single-sample enforcement for collapse is pending.

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
| `utils.nim` | Shared types: `SvType`, `Metric`, `MatchResult`, `MatchConfig`, `OutputHeader`, `formatMatchResult` |
| `intervals.nim` | `reciprocalOverlap`, `jaccard` |
| `bins.nim` | `binIndexFor`, `adjacentBins`, `TiledBuffer`, `BufferedRec` |
| `preproc.nim` | Normalize → per-(svtype,bin) temp BCF + work queue; `extraKeepInfo` (anno); `keepPassQual` (collapse); `MATCHA_BOFF` encoding |
| `matchcore.nim` | `streamJobPairs` (interval) and `streamBndJobPairs` (BND) — both return `seq[MatchPair]` (aOff, bOff, sim). Slim-BCF decode helpers. |
| `match.nim` | match-mode adapter: slim-BCF resolution → `MatchResult`, thread pools (`runMatchJobsWithPool`, `runMatchPairJobsWithPool`), TSV output |
| `anno.nim` | anno-mode: expression parser, `applyAggFunc`, slim-B DB INFO resolution, output VCF assembly |
| `mergecore.nim` | header merge, k-way slim-BCF merge, `buildSimilarityMap`, union-find, agglomerative clustering, `selectRepresentative` |
| `collapse.nim` | collapse-mode driver: per-caller preproc, two-pass matching, clustering, representative output |
| `log.nim` | Verbose logging (stderr, timestamped) |
