# sv-matcha

Compiled, efficient structural variant (SV) matching and annotation tool written in Nim using hts-nim.

## Modes

| Mode | Status |
|---|---|
| `matcha match` — pairwise matching between two SV callsets | complete |
| `matcha anno` — annotate a query callset from a population database VCF | complete |
| `matcha collapse` — cluster SVs from multiple callers within one sample, emit one representative per cluster | complete (single-sample enforcement pending) |
| `matcha merge` — merge SVs across samples into a cohort pVCF | planned |

DEL/DUP/INV match by coordinate + size (reciprocal overlap or Jaccard). BND records match by breakend proximity. INS is out of scope (silent skip). Genotypes are ignored.

## Build

```
nimble build    # → ./matcha
nimble test     # run test suite
```

Requires Nim ≥ 2.0 and hts-nim (vendored at `vendor/hts-nim/`, pinned v0.3.31).  
Fixture regeneration: `python3 tests/generate_fixtures.py` (needs `bcftools` + `bgzip`).

---

## matcha match

Find all pairwise matches between two SV callsets. Only matched pairs are emitted.

```
matcha match [options] callsetA callsetB    # cross-callset
matcha match --self [options] INPUT         # self-match (each pair once, no self-self)

Options:
  --min-overlap FLOAT    minimum reciprocal overlap (0.0–1.0)  ← exactly one required
  --min-jaccard FLOAT    minimum Jaccard index (0.0–1.0)       ←
  --bnd-slop INT         max breakend offset for BND matches (default: 100)
  --self                 match a single input against itself
  --threads INT          worker threads (default: 1)
  --tmp-dir PATH         temp directory (default: system temp)
  --output PATH          output file (default: stdout)
  -v, --verbose          verbose logging to stderr
  -h, --help
```

Inputs may be `.vcf.gz` or `.bcf`; format is auto-detected. `-v` is accepted before or after the subcommand.

### Matching logic

**DEL/DUP/INV** — candidates share CHROM and SVTYPE with sizes within the active threshold ratio.
- **Reciprocal overlap** (`--min-overlap`): `overlap / max(lenA, lenB)` — standard truvari/bedtools definition.
- **Jaccard** (`--min-jaccard`): `overlap / union`.

**BND** — candidates share CHROM and CHR2 with both breakends within `--bnd-slop`:  
`|POS_A − POS_B| < slop` and `|POS2_A − POS2_B| < slop`.  
Similarity: `(2·slop − |dPOS| − |dPOS2|) / (2·slop)`. Strand is ignored.

### Output

Tab-separated. Skip comment lines with `grep -v ^#` or `awk '!/^#/'`.

```
##matcha_metric=<overlap|jaccard>
#CHROM_A  POS_A  END_A  ID_A  CHROM_B  POS_B  END_B  ID_B  SVTYPE  SIMILARITY
```

BND rows emit `.` for `END_A` and `END_B`. `CHROM_A` and `CHROM_B` are always equal in the current matcher (each job is per-chromosome); the two columns make the schema explicit for future cross-chrom pairing.

---

## matcha anno

Annotate an input VCF/BCF with INFO fields from a database VCF, based on SV matches.

```
matcha anno [options] input database

  -a OUTFIELD=FUNC(SRCFIELD)    annotation expression (repeatable, ≥1 required)
  -o PATH                       output (.vcf | .vcf.gz | .bcf); default stdout VCF
  --min-overlap FLOAT           (exactly one of these two is required)
  --min-jaccard FLOAT
  --bnd-slop INT                default 100
  --overwrite                   allow replacing existing INFO fields
  --threads INT                 default 1
  --tmp-dir PATH
  -v, --verbose
  -h, --help
```

### Aggregation functions

`max | min | mean | first | last | best | all | unique`

- `mean` always emits Float (integer source coerces).
- `all` / `unique` emit `Number=.` lists; others emit `Number=1` scalars.
- `best` picks the match with the highest SIMILARITY (ties: earliest `posB`).
- List-valued source fields are flattened across all matches before aggregation.

### Implicit MATCHA_* fields

Available as SRCFIELD in any `-a` expression:

- `MATCHA_COUNT` — Integer, number of database matches. Unmatched records emit `0` when this is the SRCFIELD; other expressions leave OUTFIELD absent.
- `MATCHA_SIMILARITY` — Float vector, one value per match (interval metric or BND proximity score).

A `##matcha_metric=<overlap|jaccard>` line is written to the output header.

---

## matcha collapse

Cluster equivalent SVs from N single-sample callsets (e.g. Delly + Manta + GRIDSS run on the same sample) and emit one representative record per cluster, with provenance INFO fields recording which callers contributed.

```
matcha collapse [options] [Name:]callset1.bcf [Name:]callset2.bcf ...

  --min-overlap FLOAT          (exactly one of these two is required)
  --min-jaccard FLOAT
  --bnd-slop INT               default 100
  --linkage average|single|complete   agglomerative linkage (default: average)
  --priority CRITERIA          comma-separated tiebreak cascade for representative
                               selection: PASS, QUAL, CENTRE, ORDER
                               default: PASS,CENTRE,ORDER
                               ORDER is always appended as the final tiebreaker
  --format FIELDS              comma-separated FORMAT fields to carry (default: GT)
  --info FIELDS                comma-separated INFO fields to keep
                               (default: all, post conflict resolution)
  -o, --output PATH            output (.vcf | .vcf.gz | .bcf); default stdout VCF
  --threads INT                default 1
  --tmp-dir PATH
  -v, --verbose
  -h, --help
```

### Input naming

Each positional may be prefixed with `Name:` (e.g. `Delly:delly.bcf`); without a prefix, the basename without extension is used. Names appear in provenance fields and drive the `CENTRE` priority criterion.

### Representative selection

Within each cluster, the representative record is picked by walking the `--priority` cascade until one criterion decides:

- `PASS` — prefer records with `FILTER=PASS`.
- `QUAL` — prefer higher `QUAL` (BCF missing values lose).
- `CENTRE` — prefer callers earlier in the command line.
- `ORDER` — implicit final tiebreaker; preserves record order from the input.

### Linkage

`--linkage` controls the agglomerative cluster definition: `single` (any pair above threshold links), `average` (mean pairwise similarity ≥ threshold; default), or `complete` (all pairwise similarities must be ≥ threshold).

### Output

Per-cluster INFO fields added to the output:

- `SOURCE` — name of the caller providing the representative record.
- `SOURCELIST` — names of all callers contributing to the cluster (CLI order).
- `N_SOURCE` — distinct input callsets in the cluster.
- `N_MERGED` — total records merged into the cluster (≥ `N_SOURCE`).

Conflict-resolved INFO fields (e.g. differing `SVLEN` across callers) and FORMAT fields configured via `--info` / `--format` are carried through from the representative record.

---

## Dependencies

- [hts-nim](vendor/hts-nim/) — VCF/BCF I/O (vendored, pinned v0.3.31)
- `bcftools` — fixture generation only
- Nim ≥ 2.0 standard library
