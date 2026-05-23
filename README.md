# sv-matcha

Compiled, efficient structural variant (SV) matching and annotation tool written in Nim using [hts-nim](https://github.com/brentp/hts-nim) and [htslib](https://github.com/samtools/htslib).

---

## Modes

| Mode | Description | Typical usage |
|---|---|---|
| [`matcha match`](#matcha-match) | Pairwise matching between two SV callsets; emits a TSV of matched pairs | `matcha match truth.vcf.gz calls.vcf.gz > pairs.tsv` |
| [`matcha anno`](#matcha-anno) | Annotate a query callset with INFO fields from a population database VCF | `matcha anno -a AF=max(AF) calls.bcf gnomad-sv.bcf -o annotated.bcf` |
| [`matcha collapse`](#matcha-collapse) | Cluster SVs from multiple callers run on one sample; emit one representative per cluster | `matcha collapse Delly:delly.bcf Manta:manta.bcf CNVnator:cnvnator.bcf -o sample.bcf` |
| [`matcha merge`](#matcha-merge) | Merge per-sample SV callsets into a multi-sample cohort pVCF with AC/AN/AF | `matcha merge sample1.bcf sample2.bcf sample3.bcf -o cohort.bcf` |

---

## How it works

Matcha decides which structural variants from different callsets refer to the same underlying event. Two records match when they share SVTYPE and CHROM and their coordinates and sizes agree closely enough — reciprocal overlap or Jaccard for DEL/DUP/INV, breakend proximity for BND, position and size similarity for INS. Genotypes are ignored.

`match` and `anno` emit those pairwise matches directly. `collapse` and `merge` go further: they group records into clusters by agglomerative linkage, then pick one representative per cluster.

Inputs are streamed and indexed lazily, so memory stays flat in callset size and the work parallelises over chromosomes and size bins. See [DESIGN.md](DESIGN.md) for the architecture.

---

## Installation

### Precompiled binary (recommended)

Download the latest Linux x86_64 binary from the [latest release](https://github.com/jemunro/matcha/releases/latest):

```bash
curl -fsSL https://github.com/jemunro/matcha/releases/latest/download/matcha -o matcha
chmod +x matcha
mv matcha /usr/local/bin/  # or anywhere on your PATH
```

The binary dynamically links against `libhts.so` — see [Providing htslib](#providing-htslib). To build from source instead, see [Build from source](#build-from-source).

### Container image

Pre-built images bundling matcha alongside `bcftools` (which also satisfies the `libhts.so` runtime dependency) are published to GHCR for every tagged release. `:latest` tracks the most recent release.

```bash
# Docker
docker pull ghcr.io/jemunro/matcha/matcha-bcftools:latest
docker run --rm -v "$PWD:/data" -w /data ghcr.io/jemunro/matcha/matcha-bcftools \
  matcha match truth.vcf.gz calls.vcf.gz > pairs.tsv

# Apptainer / Singularity
apptainer pull docker://ghcr.io/jemunro/matcha/matcha-bcftools:latest
apptainer exec matcha-bcftools_latest.sif \
  matcha match truth.vcf.gz calls.vcf.gz > pairs.tsv
```

The image works as a drop-in container for Nextflow / Snakemake pipelines — both binaries `matcha` and `bcftools` are on `$PATH`.

---

## Quick start

```bash
# 1. Pairwise match — find SVs in calls.vcf.gz that recur in truth.vcf.gz
matcha match --min-jaccard 0.75 truth.vcf.gz calls.vcf.gz > pairs.tsv

# 2. Annotate calls with population AF from a database VCF
matcha anno -a AF=max(AF) -a AC=first(AC) calls.bcf gnomad-sv.bcf -o annotated.bcf

# 3. Collapse three caller outputs from the same sample into a unified callset
matcha collapse \
  Delly:delly.bcf Manta:manta.bcf CNVnator:cnvnator.bcf \
  --priority PASS,QUAL,CENTRE,ORDER -o sample.bcf

# 4. Merge per-sample callsets into a multi-sample cohort pVCF
matcha merge sample1.bcf sample2.bcf sample3.bcf --missing-to-ref -o cohort.bcf
```

---

## Matching semantics

All four subcommands share the same matching engine and the same threshold options. SVTYPE-specific rules:

**DEL / DUP / INV** — candidates share CHROM and SVTYPE with sizes within the active threshold ratio.
- **Reciprocal overlap** (`--min-overlap`): `overlap / max(lenA, lenB)` — the standard truvari/bedtools definition.
- **Jaccard** (`--min-jaccard`): `overlap / union`.

The two metrics are mutually exclusive; the default is `--min-jaccard 0.75`.

**BND** — candidates share CHROM and CHR2 with both breakends within `--bnd-slop` (default 50):
`|POS_A − POS_B| < slop` and `|POS2_A − POS2_B| < slop`. Similarity is `(2·slop − |dPOS| − |dPOS2|) / (2·slop)`. Strand is ignored.

**INS** — candidates share CHROM with `|POS_A − POS_B| < --ins-slop` (default 50). Similarity combines position and size:
- `pos_sim = (slop − |dPOS|) / slop`
- `len_sim = min(SVLEN_A, SVLEN_B) / max(SVLEN_A, SVLEN_B)`
- `sim = sqrt(pos_sim · len_sim)` — must be ≥ `--min-ins-sim` (default 0.75).

SVLEN for INS is resolved from the first available of: `INFO/INSLEN`, `INFO/SVLEN`, `len(ALT) − len(REF)` (sequence-resolved ALTs only), `len(LEFT_SVINSSEQ) + len(RIGHT_SVINSSEQ)`. Records with no resolvable length are skipped with reason `unresolvable_ins_len`.

`TRA` records are not supported; they are warned and skipped.

---

## matcha match

Find all pairwise matches between two SV callsets. Only matched pairs are emitted.

```
matcha match [options] callsetA callsetB    # cross-callset
matcha match --self [options] INPUT         # self-match (each pair once, no self-self)

Options:
  --min-overlap FLOAT    minimum reciprocal overlap (0.0–1.0)  ← mutually exclusive;
  --min-jaccard FLOAT    minimum Jaccard index (0.0–1.0)       ← default --min-jaccard 0.75
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

Inputs may be `.vcf.gz` or `.bcf`; format is auto-detected. `-v` is accepted before or after the subcommand. See [Matching semantics](#matching-semantics) for how candidates are evaluated.

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
  --min-overlap FLOAT           (mutually exclusive; default: --min-jaccard 0.75)
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

See [Matching semantics](#matching-semantics) for matching rules.

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

  --min-overlap FLOAT          (mutually exclusive; default: --min-jaccard 0.75)
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
                               (default: SVTYPE,SVLEN,END,CHR2,POS2 only)
  -o, --output PATH            output (.vcf | .vcf.gz | .bcf); default stdout VCF
  --chrs CHR[,CHR...]          restrict to listed chromosomes (filters records + header contigs)
  --threads INT                default 1
  --tmp-dir PATH
  -v, --verbose
  -h, --help
```

See [Matching semantics](#matching-semantics) for how cluster candidates are evaluated.

### Input naming

Each positional may be prefixed with `Name:` (e.g. `Delly:delly.bcf`); without a prefix, the basename without extension is used. Names appear in provenance fields and drive the `ORDER` priority criterion.

### Representative selection

Within each cluster, the representative record is picked by walking the `--priority` cascade until one criterion decides:

- `PASS` — prefer records with `FILTER=PASS`.
- `QUAL` — prefer higher `QUAL` (BCF missing values lose).
- `CENTRE` — prefer records closest to the cluster centre, i.e. those with the highest mean pairwise similarity to other cluster members.
- `ORDER` — prefer records from callers listed earlier on the command line (via `CALLER_IDX`). Always appended as the final criterion in the cascade.

Any remaining ties are broken deterministically by taking the first record in input order (lowest `SRC_INDEX`, which for records from a single caller corresponds to the lowest position since inputs are coordinate-sorted).

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
matcha merge [options] callset1.bcf callset2.bcf ...

  --min-overlap FLOAT           (mutually exclusive; default: --min-jaccard 0.75)
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

See [Matching semantics](#matching-semantics) for how cluster candidates are evaluated.

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

## Build from source

```
git clone --recurse-submodules https://github.com/jemunro/matcha
cd matcha
nimble build    # → ./matcha
nimble release  # optimised, stripped
nimble test     # run test suite
```

Requirements: Nim ≥ 2.0, a C compiler, and hts-nim (vendored at `vendor/hts-nim/`, pinned v0.3.31), which links against htslib ≥ 1.10. Tests additionally require `bcftools` and `bgzip` for fixture generation (`python3 tests/generate_fixtures.py`).

## Providing htslib

`hts-nim` links against the system `libhts.so` at runtime. Pick whichever route fits your environment:

- **conda/mamba**: `mamba install -c bioconda 'htslib>=1.10'` — gets you a recent build (`libhts.so` lands in `$CONDA_PREFIX/lib`).
- **Debian/Ubuntu**: `apt install libhts-dev` — check the version with `apt-cache policy libhts-dev`; Ubuntu ≥ 22.04 ships ≥ 1.13. Older releases need a backport or source build.
- **macOS**: `brew install htslib`.
- **From source**: clone `https://github.com/samtools/htslib`, check out a `1.x` tag (≥ 1.10), then `autoreconf -i && ./configure && make && sudo make install`.

If `libhts.so` is not on the default loader path, point to it at runtime with `LD_LIBRARY_PATH=/path/to/htslib/lib ./matcha …`. Confirm which library is being picked up with `ldd ./matcha | grep hts` (Linux) or `otool -L ./matcha | grep hts` (macOS).

## Citations

> Brent S Pedersen, Aaron R Quinlan. **hts-nim: scripting high-performance genomic analyses.** *Bioinformatics*, Volume 34, Issue 19, October 2018, Pages 3387–3389. <https://doi.org/10.1093/bioinformatics/bty358>
