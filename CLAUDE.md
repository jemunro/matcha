# matcha

Compiled SV matching and annotation tool (Nim + hts-nim). See [README.md](README.md) for user docs, [DESIGN.md](DESIGN.md) for architecture detail.

## Status

Milestones 1–2 (`matcha match`, `matcha anno`) complete. SVTYPE=INS support is out of scope (silent skip). `matcha collapse` and `matcha merge` are planned.

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
| `utils.nim` | Shared types: `SvType`, `Metric`, `MatchResult`, `MatchConfig` |
| `intervals.nim` | `reciprocalOverlap`, `jaccard` |
| `bins.nim` | `binIndexFor`, `adjacentBins`, `TiledBuffer` |
| `preproc.nim` | Normalize → per-(svtype,bin) temp BCF + work queue; `extraKeepInfo` for anno |
| `matchcore.nim` | `streamJobPairs[B,R]` (interval) and `streamBndJobPairs[B,R]` (BND) |
| `match.nim` | match-mode adapter: thread pool, self-mode dedup, TSV output |
| `anno.nim` | anno-mode: expression parser, `applyAggFunc`, output VCF assembly |
| `log.nim` | Verbose logging (stderr, timestamped) |
