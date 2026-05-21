# sv-matcha

Compiled, efficient structural variant (SV) matching and annotation tool written in Nim using hts-nim.

## Modes

| Mode | Status |
|---|---|
| `matcha match` — pairwise matching between two SV callsets | complete |
| `matcha anno` — annotate a query callset from a population database VCF | complete |
| `matcha collapse` — cluster SVs from multiple callers within one sample, emit one representative per cluster | complete |
| `matcha merge` — merge SVs across samples into a cohort pVCF | complete |

DEL/DUP/INV match by coordinate + size (reciprocal overlap or Jaccard). BND records match by breakend proximity. INS records match by position proximity plus size ratio. Genotypes are ignored.

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
  --bnd-slop INT         max breakend offset for BND matches (default: 50)
  --min-ins-sim FLOAT    minimum INS combined sim = sqrt(pos_sim·len_sim) (default: 0.75)
  --ins-slop INT         max position offset for INS matches (default: 50)
  --self                 match a single input against itself
  --info FIELDS          comma-separated INFO fields to include as INFO_A/INFO_B columns
  --chrs CHR[,CHR...]    restrict to listed chromosomes (filters records; no header change for TSV output)
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

**INS** — candidates share CHROM with `|POS_A − POS_B| < --ins-slop`. Similarity combines position and size:
- `pos_sim = (slop − |dPOS|) / slop`
- `len_sim = min(SVLEN_A, SVLEN_B) / max(SVLEN_A, SVLEN_B)`
- `sim = sqrt(pos_sim · len_sim)` — must be ≥ `--min-ins-sim`.

SVLEN is resolved from the first available of: `INFO/INSLEN`, `INFO/SVLEN`, `len(ALT) − len(REF)` (sequence-resolved ALTs only), `len(LEFT_SVINSSEQ) + len(RIGHT_SVINSSEQ)`. Records with no resolvable length are skipped with reason `unresolvable_ins_len`.

### Output

Tab-separated. Skip comment lines with `grep -v ^#` or `awk '!/^#/'`.

```
##matcha_metric=<overlap|jaccard>
#CHROM_A  POS_A  ID_A  CHROM_B  POS_B  ID_B  SVTYPE  SIMILARITY
```

`CHROM_A` and `CHROM_B` are always equal (each job is per-chromosome).

With `--info SVLEN,END,AF` two extra columns are inserted — `INFO_A` after `ID_A` and `INFO_B` after `ID_B` — containing the requested fields in VCF INFO format (`KEY=VALUE;KEY=VALUE`). Fields absent on a record are omitted from that cell; `.` is emitted when none are present.

```
#CHROM_A  POS_A  ID_A  INFO_A  CHROM_B  POS_B  ID_B  INFO_B  SVTYPE  SIMILARITY
chr1      1000   DEL_A_01  SVLEN=-1000;END=2000  chr1  1000  DEL_B_01  SVLEN=-1000;END=2000  DEL  1.000000
```

---

## matcha anno

Annotate an input VCF/BCF with INFO fields from a database VCF, based on SV matches.

```
matcha anno [options] input database

  -a OUTFIELD=FUNC(SRCFIELD)    annotation expression (repeatable, ≥1 required)
  -o PATH                       output (.vcf | .vcf.gz | .bcf); default stdout VCF
  --min-overlap FLOAT           (exactly one of these two is required)
  --min-jaccard FLOAT
  --bnd-slop INT                default 50
  --min-ins-sim FLOAT           default 0.75
  --ins-slop INT                default 50
  --overwrite                   allow replacing existing INFO fields
  --chrs CHR[,CHR...]           restrict to listed chromosomes (filters records + header contigs)
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
  --bnd-slop INT               default 50
  --min-ins-sim FLOAT          default 0.75
  --ins-slop INT               default 50
  --linkage average|single|complete   agglomerative linkage (default: average)
  --priority CRITERIA          comma-separated tiebreak cascade for representative
                               selection: PASS, QUAL, CENTRE, ORDER
                               default: PASS,CENTRE,ORDER
                               ORDER is always appended as the final tiebreaker
  --format FIELDS              comma-separated FORMAT fields to carry (default: GT)
  --info FIELDS                comma-separated INFO fields to keep
                               (default: all, post conflict resolution)
  -o, --output PATH            output (.vcf | .vcf.gz | .bcf); default stdout VCF
  --chrs CHR[,CHR...]          restrict to listed chromosomes (filters records + header contigs)
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

- `CALLERS` — caller names in the cluster: representative first, followed by others in CLI order.
- `N_CALLERS` — distinct input callsets in the cluster.
- `N_MERGED` — total records merged into the cluster (≥ `N_CALLERS`).

Conflict-resolved INFO fields (e.g. differing `SVLEN` across callers) and FORMAT fields configured via `--info` / `--format` are carried through from the representative record. Conflict rules:

- **Compatible** (same Number + Type across all callers that declare the field): single merged def, no rename.
- **Number-only conflict** (same Type, different Number — e.g. `Number=1` vs `Number=.`): silently widened to `Number=.` in the output. No rename, no warning.
- **Type conflict** (e.g. one caller declares `Integer`, another `String`): each caller's instance is renamed to `FIELD_<callerName>` in the output. A warning is logged only when the field is in the `--info` / `--format` filter; warnings for fields the user didn't request are suppressed.

---

## matcha merge

Merge per-sample SV callsets (typically `matcha collapse` outputs) into a single
multi-sample cohort pVCF: one row per cluster, per-sample FORMAT columns,
cohort INFO (AC/AN/AF) computed from the assembled GTs.

```
matcha merge [options] [Name:]callset1.bcf [Name:]callset2.bcf ...

  --min-overlap FLOAT           (exactly one of these two is required)
  --min-jaccard FLOAT
  --bnd-slop INT                default 50
  --min-ins-sim FLOAT           default 0.75
  --ins-slop INT                default 50
  --linkage average|single|complete   agglomerative linkage (default: average)
  --priority CRITERIA           tiebreak cascade for representative selection
                                default: PASS,CENTRE,ORDER
  --format FIELDS               comma-separated FORMAT fields to carry per sample
                                default: GT (auto-added if absent)
  --info FIELDS                 comma-separated INFO fields to keep from representative
                                default: only auto-extracted + cohort + CALLERS
  -o, --output PATH             output (.vcf | .vcf.gz | .bcf); default stdout VCF
  --chrs CHR[,CHR...]           restrict to listed chromosomes (filters records + header contigs)
  --missing-to-ref              treat absent samples as 0/0 (count toward AN; like bcftools merge)
  --threads INT                 default 1
  --tmp-dir PATH
  -v, --verbose
  -h, --help
```

### Input requirements

- ≥ 2 input files.
- Each input must have **exactly 1 sample column**; multi-sample inputs are rejected.
- All sample IDs across inputs must be **distinct**; collisions are rejected.

### Output

- N sample columns named with each input's own sample ID, in CLI order.
- Per-cluster cohort INFO fields:
  - `AC` (`Number=A,Integer`) — alt allele count across called genotypes.
  - `AN` (`Number=1,Integer`) — total alleles called (sums over non-missing GTs).
  - `AF` (`Number=A,Float`) — `AC / AN`; missing (`AF=.`) when `AN == 0`.
- When any input record carries `INFO/CALLERS`, the output also emits:
  - `CALLERS` — union across cluster members (representative caller first).
  - `N_CALLERS` — distinct count.
- Missing samples are written as `GT=./.` with missing values for other carried FORMAT fields. With `--missing-to-ref`, absent samples are instead written as `GT=0/0` and counted toward `AN` (mirroring `bcftools merge --missing-to-ref`); the AF denominator changes accordingly. The flag affects absent samples only — present samples with in-call missing alleles (e.g. `./1`) are untouched.
- For BND records, the original bracket-form ALT (`N[chr:pos[` etc.) is preserved.
- For INS records, the original sequence-resolved ALT (e.g. from Delly or sequence-resolved Manta calls) is preserved on the representative record when present; symbolic `<INS>` is used otherwise.

### Notes

- Same-sample collisions inside a single cluster are **not blocked** — if two records from
  the same sample happen to cluster together, the priority cascade (PASS → QUAL → ORDER)
  picks one per sample column and a throttled warning is emitted.
- `GT` is silently added to `--format` when absent so cohort AC/AN/AF can be computed.
- `N_MERGED` (collapse-specific) is not emitted.

---

## Dependencies

- [hts-nim](vendor/hts-nim/) — VCF/BCF I/O (vendored, pinned v0.3.31)
- `bcftools` — fixture generation only
- Nim ≥ 2.0 standard library
