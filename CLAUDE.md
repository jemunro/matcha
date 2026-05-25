# matcha

Compiled SV matching and annotation tool (Nim + hts-nim). See [README.md](README.md) for user docs, [DESIGN.md](DESIGN.md) for architecture detail.

## Scope

All four subcommands (`match`, `anno`, `collapse`, `merge`) support DEL/DUP/INV/BND/INS. TRA is warned and skipped.

## Build & test

```
nimble build          # → ./matcha
nimble test           # run test suite
python3 tests/generate_fixtures.py   # regenerate fixtures (needs bcftools + bgzip)
```

## Code layout

`main.nim` parses CLI args and dispatches to one of four per-mode drivers (`match.nim`, `anno.nim`, `collapse.nim`, `merge.nim`). All four call into shared `preproc.nim` (normalize + slim each input into per-(svtype, bin) temp BCFs) and `matchcore.nim` (three streamers — interval / BND / INS — each returning `seq[MatchPair]`). `mergecore.nim` provides header merging and agglomerative clustering shared by `collapse` and `merge`. See [DESIGN.md](DESIGN.md) for the full module map and architecture notes.
