#!/usr/bin/env python3
"""Generate test fixtures for matcha.

Usage: python tests/generate_fixtures.py [--output-dir tests/fixtures]
Requires: bgzip, bcftools on PATH.
"""
import argparse
import os
import subprocess
import sys

HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chrX,length=156040895>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

# DB fixture uses an augmented header carrying annotation INFO fields.
# Match tests don't touch this; anno tests pull values off these records.
DB_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency in the database cohort">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count in the database cohort">
##INFO=<ID=CALLER,Number=1,Type=String,Description="Single primary caller">
##INFO=<ID=CALLERS,Number=.,Type=String,Description="All callers reporting this SV">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chrX,length=156040895>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

# Records are dicts so edge cases (missing fields, conflicting fields, missing IDs)
# can be expressed without contorting a fixed tuple shape. Each dict carries:
#   - VCF columns: chrom, pos, id, ref, alt, info (preformatted INFO string)
#   - Matching truth: svtype, svlen, end (None if record is intentionally
#     malformed or unmatchable; compute_expected then skips it)
#   - For BND records: chr2, pos2 are added so compute_bnd_expected can pair them.
def normal(chrom, pos, rid, ref, alt, svtype, svlen, end):
    return {
        "chrom": chrom, "pos": pos, "id": rid, "ref": ref, "alt": alt,
        "info": f"SVTYPE={svtype};SVLEN={svlen};END={end}",
        "svtype": svtype, "svlen": svlen, "end": end,
    }


def bnd(chrom, pos, rid, chr2, pos2):
    """BND record using the t[p[ ALT form. matcha parses the bracket; INFO
    only carries SVTYPE=BND. svlen=0 / end=pos+1 are sentinels matching
    matcha's slim-BCF encoding (not used for matching)."""
    alt = f"N[{chr2}:{pos2}["
    return {
        "chrom": chrom, "pos": pos, "id": rid, "ref": "N", "alt": alt,
        "info": "SVTYPE=BND",
        "svtype": "BND", "svlen": 0, "end": pos + 1,
        "chr2": chr2, "pos2": pos2,
    }


RECORDS_A = [
    # TC01 exact match
    normal("chr1", 1000,  "DEL_A_01", "N", "<DEL>", "DEL",  -1000, 2000),
    # TC02 partial overlap above threshold (shift +100)
    normal("chr1", 3000,  "DEL_A_02", "N", "<DEL>", "DEL",  -1000, 4000),
    # TC03 partial overlap below threshold (shift +600 in B)
    normal("chr1", 5000,  "DEL_A_03", "N", "<DEL>", "DEL",  -1000, 6000),
    # TC04 no overlap
    normal("chr1", 7000,  "DEL_A_04", "N", "<DEL>", "DEL",  -1000, 8000),
    # TC05 SVTYPE mismatch (DEL in A, DUP in B)
    normal("chr1", 11000, "DEL_A_05", "N", "<DEL>", "DEL",  -1000, 12000),
    # TC06 multiple matches for one A record
    normal("chr1", 13000, "DEL_A_06", "N", "<DEL>", "DEL",  -1000, 14000),
    # TC07 unmatched A (no B record nearby)
    normal("chr1", 15000, "DEL_A_07", "N", "<DEL>", "DEL",  -1000, 16000),
    # TC08 size asymmetry (large A, small B)
    normal("chr1", 17000, "DEL_A_08", "N", "<DEL>", "DEL",  -5000, 22000),
    # TC09 multi-chromosome
    normal("chr2", 1000,  "DUP_A_09", "N", "<DUP>", "DUP",   2000, 3000),
    normal("chrX", 5000,  "INV_A_09", "N", "<INV>", "INV",   2000, 7000),
    # TC10 all three SVTYPEs on chr1
    normal("chr1", 25000, "DUP_A_10", "N", "<DUP>", "DUP",   2000, 27000),
    normal("chr1", 28000, "INV_A_10", "N", "<INV>", "INV",   2000, 30000),
    # TC11 BND (matchable) + INS (silently skipped)
    bnd("chr1", 31000, "BND_A_11", chr2="chr1", pos2=32000),
    {"chrom": "chr1", "pos": 32000, "id": "INS_A_11", "ref": "N",
     "alt": "<INS>", "info": "SVTYPE=INS;SVLEN=100;END=32100",
     "svtype": None, "svlen": None, "end": None},
    # TC12 SVTYPE only in ALT (no INFO/SVTYPE) → kept, classified as DEL
    {"chrom": "chr1", "pos": 33000, "id": "DEL_A_12_alt_only", "ref": "N",
     "alt": "<DEL>", "info": "SVLEN=-1000;END=34000",
     "svtype": "DEL", "svlen": -1000, "end": 34000},
    # TC13 SVTYPE conflict (INFO=DUP, ALT=<DEL>) → ALT wins, kept as DEL
    {"chrom": "chr1", "pos": 35000, "id": "DEL_A_13_conflict", "ref": "N",
     "alt": "<DEL>", "info": "SVTYPE=DUP;SVLEN=-1000;END=36000",
     "svtype": "DEL", "svlen": -1000, "end": 36000},
    # TC14 missing ID → synthesized as CHROM_POS_SVTYPE_LINENUMBER
    {"chrom": "chr1", "pos": 37000, "id": ".", "ref": "N",
     "alt": "<DEL>", "info": "SVTYPE=DEL;SVLEN=-1000;END=38000",
     "svtype": "DEL", "svlen": -1000, "end": 38000},
    # TC15 missing END and SVLEN → skipped (skUnresolvableEnd)
    {"chrom": "chr1", "pos": 39000, "id": "DEL_A_15_no_end_no_svlen", "ref": "N",
     "alt": "<DEL>", "info": "SVTYPE=DEL",
     "svtype": None, "svlen": None, "end": None},
    # TC16 END < POS → skipped (skEndLePos)
    {"chrom": "chr1", "pos": 41000, "id": "DEL_A_16_bad_end", "ref": "N",
     "alt": "<DEL>", "info": "SVTYPE=DEL;SVLEN=-500;END=40500",
     "svtype": None, "svlen": None, "end": None},
    # TC17 END/SVLEN inconsistent (>10%) → kept, normalized to SVLEN=END-POS=1000
    {"chrom": "chr1", "pos": 43000, "id": "DEL_A_17_inconsistent", "ref": "N",
     "alt": "<DEL>", "info": "SVTYPE=DEL;SVLEN=-1500;END=44000",
     "svtype": "DEL", "svlen": -1000, "end": 44000},
    # TC18 — self-mode cluster: three DELs that overlap each other (chr1).
    # Pairwise reciprocal overlaps at threshold 0.5 (all kept):
    #   (18a, 18b): 900/1000 = 0.90
    #   (18a, 18c): 950/1000 = 0.95
    #   (18b, 18c): 950/1000 = 0.95
    # No B-fixture counterparts → invisible in cross-callset tests.
    normal("chr1", 50000, "DEL_A_18a", "N", "<DEL>", "DEL", -1000, 51000),
    normal("chr1", 50100, "DEL_A_18b", "N", "<DEL>", "DEL", -1000, 51100),
    normal("chr1", 50050, "DEL_A_18c", "N", "<DEL>", "DEL", -1000, 51050),
    # TC19 — self-mode pair on chr2 (different chrom, different svtype).
    # (19a, 19b): 1500/2000 = 0.75.
    normal("chr2", 5000, "DUP_A_19a", "N", "<DUP>", "DUP", 2000, 7000),
    normal("chr2", 5500, "DUP_A_19b", "N", "<DUP>", "DUP", 2000, 7500),
    # TC20 — inter-chrom BND identical mate. With B identical, sim=1.0.
    bnd("chr1", 60000, "BND_A_20", chr2="chr2", pos2=70000),
    # TC21 — slight offset on both breakends. dPOS=50, dPOS2=30.
    # sim at slop=100 = (200-50-30)/200 = 0.60.  At slop=20: |50|>=20 → reject.
    bnd("chr1", 60500, "BND_A_21", chr2="chr2", pos2=70000),
    # TC22 — no B counterpart on chr1:80000 → unmatched.
    bnd("chr1", 80000, "BND_A_22", chr2="chr2", pos2=90000),
    # TC23 — CHR2 mismatch: A has chrX mate; B has chr2 mate. → no match.
    bnd("chr1", 62000, "BND_A_23", chr2="chrX", pos2=70000),
    # TC24 — malformed BND ALT (truncated). → warn-skip (skMalformedBnd).
    {"chrom": "chr1", "pos": 64000, "id": "BND_A_24", "ref": "N",
     "alt": "N[chr1:", "info": "SVTYPE=BND",
     "svtype": None, "svlen": None, "end": None},
    # TC25 — TRA record. → warn-skip (skUnsupportedTra).
    # No CHR2 in INFO since the fixture header doesn't declare it; bcftools
    # rejects undeclared INFO fields. matcha warn-skips TRA regardless of
    # mate-position presence.
    {"chrom": "chr1", "pos": 70000, "id": "TRA_A_25", "ref": "N",
     "alt": "<TRA>", "info": "SVTYPE=TRA;END=70001",
     "svtype": None, "svlen": None, "end": None},
    # TC27 — intra-chrom BND mate pair on chr2 (for --self mode).
    # 27a/27b: dPOS=20, dPOS2=20 → sim at slop=100 = (200-20-20)/200 = 0.80.
    bnd("chr2", 60000, "BND_A_27a", chr2="chr2", pos2=80000),
    bnd("chr2", 60020, "BND_A_27b", chr2="chr2", pos2=80020),
]

RECORDS_B = [
    # TC01 exact match
    normal("chr1", 1000,  "DEL_B_01", "N", "<DEL>", "DEL",  -1000, 2000),
    # TC02 partial overlap (shifted +100)
    normal("chr1", 3100,  "DEL_B_02", "N", "<DEL>", "DEL",  -1000, 4100),
    # TC03 below threshold (shifted +600)
    normal("chr1", 5600,  "DEL_B_03", "N", "<DEL>", "DEL",  -1000, 6600),
    # TC04 no overlap (far away)
    normal("chr1", 9000,  "DEL_B_04", "N", "<DEL>", "DEL",  -1000, 10000),
    # TC05 SVTYPE mismatch (DUP in B)
    normal("chr1", 11000, "DUP_B_05", "N", "<DUP>", "DUP",   1000, 12000),
    # TC06 two B records that match DEL_A_06
    normal("chr1", 13050, "DEL_B_06a","N", "<DEL>", "DEL",  -1000, 14050),
    normal("chr1", 12950, "DEL_B_06b","N", "<DEL>", "DEL",  -1000, 13950),
    # TC07: no B record for DEL_A_07 (intentionally absent)
    # TC08 small B for large A
    normal("chr1", 17000, "DEL_B_08", "N", "<DEL>", "DEL",  -1000, 18000),
    # TC09 multi-chromosome
    normal("chr2", 1000,  "DUP_B_09", "N", "<DUP>", "DUP",   2000, 3000),
    normal("chrX", 5000,  "INV_B_09", "N", "<INV>", "INV",   2000, 7000),
    # TC10 all three SVTYPEs
    normal("chr1", 25000, "DUP_B_10", "N", "<DUP>", "DUP",   2000, 27000),
    normal("chr1", 28000, "INV_B_10", "N", "<INV>", "INV",   2000, 30000),
    # TC11 BND in B — identical mate. Matches BND_A_11 at sim=1.0.
    bnd("chr1", 31000, "BND_B_11", chr2="chr1", pos2=32000),
    # TC20 inter-chrom BND, identical to A → sim=1.0
    bnd("chr1", 60000, "BND_B_20", chr2="chr2", pos2=70000),
    # TC21 inter-chrom BND, offset by 50/30 → sim=0.60 at default slop
    bnd("chr1", 60550, "BND_B_21", chr2="chr2", pos2=70030),
    # TC23 — A's BND_A_23 has chrX mate; this B record has chr2 mate.
    # Position matches (chr1:62000) but CHR2 mismatch → must NOT pair.
    bnd("chr1", 62000, "BND_B_23", chr2="chr2", pos2=70000),
]


def recip_overlap(pA, eA, pB, eB):
    ovl = min(eA, eB) - max(pA, pB)
    if ovl <= 0:
        return 0.0
    return ovl / max(eA - pA, eB - pB)


def jaccard(pA, eA, pB, eB):
    ovl = min(eA, eB) - max(pA, pB)
    if ovl <= 0:
        return 0.0
    union = max(eA, eB) - min(pA, pB)
    return ovl / union


def db(chrom, pos, rid, alt, svtype, svlen, end, af, ac, caller, callers):
    """Build a DB-style record with anno INFO fields baked into the info string."""
    info = (
        f"SVTYPE={svtype};SVLEN={svlen};END={end};"
        f"AF={af};AC={ac};CALLER={caller};CALLERS={callers}"
    )
    return {
        "chrom": chrom, "pos": pos, "id": rid, "ref": "N", "alt": alt,
        "info": info,
        "svtype": svtype, "svlen": svlen, "end": end,
    }


def bnd_db(chrom, pos, rid, chr2, pos2, af, ac, caller, callers):
    """BND-style DB record with anno INFO fields."""
    alt = f"N[{chr2}:{pos2}["
    info = (
        f"SVTYPE=BND;"
        f"AF={af};AC={ac};CALLER={caller};CALLERS={callers}"
    )
    return {
        "chrom": chrom, "pos": pos, "id": rid, "ref": "N", "alt": alt,
        "info": info,
        "svtype": "BND", "svlen": 0, "end": pos + 1,
        "chr2": chr2, "pos2": pos2,
    }


# DB fixture records, each carrying AF/AC/CALLER/CALLERS. Designed so that
# specific A records have predictable annotation outputs:
#   DEL_A_01 (chr1:1000-2000) -> one match (DEL_DB_01)
#   DEL_A_02 (chr1:3000-4000) -> one match (DEL_DB_02)  (partial overlap)
#   DEL_A_06 (chr1:13000-14000) -> two matches (DEL_DB_06a, DEL_DB_06b)
#   DEL_A_07 (chr1:15000-16000) -> zero matches (no nearby DB record)
#   BND_A_11 (chr1:31000 → chr1:32000) -> one match (BND_DB_11)
RECORDS_DB = [
    db("chr1", 1000,  "DEL_DB_01",  "<DEL>", "DEL", -1000, 2000,
       0.10, 5,  "manta", "manta,delly"),
    db("chr1", 3000,  "DEL_DB_02",  "<DEL>", "DEL", -1000, 4000,
       0.50, 25, "delly", "delly"),
    # Two matches for DEL_A_06; B 06b is at smaller POS, so it sorts first
    # in posB order — first()/last() should be deterministic on that ordering.
    db("chr1", 12950, "DEL_DB_06b", "<DEL>", "DEL", -1000, 13950,
       0.30, 15, "delly", "delly,gatk"),
    db("chr1", 13050, "DEL_DB_06a", "<DEL>", "DEL", -1000, 14050,
       0.20, 10, "manta", "manta"),
    # BND match for BND_A_11.
    bnd_db("chr1", 31000, "BND_DB_11", chr2="chr1", pos2=32000,
           af=0.40, ac=20, caller="manta", callers="manta,gridss"),
]


def records_to_vcf(records, header=HEADER):
    lines = [header]
    for r in records:
        qual   = r.get("qual", ".")
        filt   = r.get("filter", "PASS")
        lines.append(
            f"{r['chrom']}\t{r['pos']}\t{r['id']}\t{r['ref']}\t{r['alt']}"
            f"\t{qual}\t{filt}\t{r['info']}"
        )
    return "\n".join(lines) + "\n"


def write_vcf_gz(records, path, header=HEADER):
    vcf_text = records_to_vcf(records, header)
    # bcftools sort bgzips and sorts in one step
    subprocess.run(
        ["bcftools", "sort", "-O", "z", "-o", path, "-"],
        input=vcf_text.encode(),
        check=True,
        capture_output=True,
    )
    # TBI index (used by hts-nim for VCF.gz region queries)
    subprocess.run(["bcftools", "index", path], check=True, capture_output=True)


def write_bcf(vcf_gz_path, bcf_path):
    """Convert an existing sorted vcf.gz to bcf and CSI-index it."""
    subprocess.run(
        ["bcftools", "view", "-O", "b", "-o", bcf_path, vcf_gz_path],
        check=True,
        capture_output=True,
    )
    subprocess.run(["bcftools", "index", bcf_path], check=True, capture_output=True)


def interval_matchable(r):
    return (r["svtype"] is not None and r["svtype"] != "BND"
            and r["svlen"] is not None and r["end"] is not None)


def bnd_matchable(r):
    return r["svtype"] == "BND" and r.get("chr2") and r.get("pos2") is not None


def compute_expected(records_a, records_b, metric, threshold):
    """Interval-mode expected rows. metric is 'overlap' or 'jaccard'.

    Records that aren't interval-matchable (None svtype, BND, malformed)
    are skipped — BND is handled separately by compute_bnd_expected.
    """
    assert metric in ("overlap", "jaccard")

    b_by_key = {}
    for rec in records_b:
        if not interval_matchable(rec):
            continue
        b_by_key.setdefault((rec["chrom"], rec["svtype"]), []).append(rec)

    rows = []
    for rec_a in records_a:
        if not interval_matchable(rec_a):
            continue
        chrom, posA, idA, svtype, endA = (
            rec_a["chrom"], rec_a["pos"], rec_a["id"], rec_a["svtype"], rec_a["end"],
        )
        for rec_b in b_by_key.get((chrom, svtype), []):
            posB, idB, endB = rec_b["pos"], rec_b["id"], rec_b["end"]
            sim = (recip_overlap(posA, endA, posB, endB) if metric == "overlap"
                   else jaccard(posA, endA, posB, endB))
            if sim < threshold:
                continue
            rows.append([
                chrom, str(posA), str(endA), idA,
                str(posB), str(endB), idB,
                svtype, f"{sim:.6f}",
            ])
    return rows


def compute_self_expected(records, metric, threshold):
    """Interval self-mode expected rows. metric is 'overlap' or 'jaccard'."""
    assert metric in ("overlap", "jaccard")

    by_key = {}
    for rec in records:
        if not interval_matchable(rec):
            continue
        by_key.setdefault((rec["chrom"], rec["svtype"]), []).append(rec)

    rows = []
    for (chrom, svtype), recs in by_key.items():
        recs_sorted = sorted(recs, key=lambda r: r["pos"])
        for i in range(len(recs_sorted)):
            for j in range(i + 1, len(recs_sorted)):
                a, b = recs_sorted[i], recs_sorted[j]
                sim = (recip_overlap(a["pos"], a["end"], b["pos"], b["end"])
                       if metric == "overlap"
                       else jaccard(a["pos"], a["end"], b["pos"], b["end"]))
                if sim < threshold:
                    continue
                rows.append([
                    chrom, str(a["pos"]), str(a["end"]), a["id"],
                    str(b["pos"]), str(b["end"]), b["id"],
                    svtype, f"{sim:.6f}",
                ])
    return rows


def compute_bnd_expected(records_a, records_b, slop=100):
    """BND expected rows. Two BNDs match iff same CHROM_A, same CHR2, and
    both breakends are within slop (strict <). Similarity is
    (2*slop - |dPOS| - |dPOS2|) / (2*slop)."""
    b_by_chrom = {}
    for rec in records_b:
        if not bnd_matchable(rec):
            continue
        b_by_chrom.setdefault(rec["chrom"], []).append(rec)

    rows = []
    twoSlop = 2 * slop
    for ra in records_a:
        if not bnd_matchable(ra):
            continue
        for rb in b_by_chrom.get(ra["chrom"], []):
            d1 = abs(ra["pos"] - rb["pos"])
            if d1 >= slop:
                continue
            if ra["chr2"] != rb["chr2"]:
                continue
            d2 = abs(ra["pos2"] - rb["pos2"])
            if d2 >= slop:
                continue
            sim = (twoSlop - d1 - d2) / twoSlop
            if sim <= 0:
                continue
            rows.append([
                ra["chrom"], str(ra["pos"]), ".", ra["id"],
                str(rb["pos"]), ".", rb["id"], "BND",
                f"{sim:.6f}",
            ])
    return rows


def compute_self_bnd_expected(records, slop=100):
    by_chrom = {}
    for rec in records:
        if not bnd_matchable(rec):
            continue
        by_chrom.setdefault(rec["chrom"], []).append(rec)

    rows = []
    twoSlop = 2 * slop
    for chrom, recs in by_chrom.items():
        recs_sorted = sorted(recs, key=lambda r: r["pos"])
        for i in range(len(recs_sorted)):
            for j in range(i + 1, len(recs_sorted)):
                a, b = recs_sorted[i], recs_sorted[j]
                d1 = abs(a["pos"] - b["pos"])
                if d1 >= slop:
                    continue
                if a["chr2"] != b["chr2"]:
                    continue
                d2 = abs(a["pos2"] - b["pos2"])
                if d2 >= slop:
                    continue
                sim = (twoSlop - d1 - d2) / twoSlop
                if sim <= 0:
                    continue
                rows.append([
                    chrom, str(a["pos"]), ".", a["id"],
                    str(b["pos"]), ".", b["id"], "BND",
                    f"{sim:.6f}",
                ])
    return rows


def write_expected_tsv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Collapse fixtures
# ---------------------------------------------------------------------------

COLLAPSE_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality call">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chrX,length=156040895>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

def col(chrom, pos, rid, svtype, svlen, end, qual=".", filt="PASS"):
    """Collapse fixture record."""
    sign = -1 if svtype == "DEL" else 1
    alt  = f"<{svtype}>"
    info = f"SVTYPE={svtype};SVLEN={sign * abs(svlen)};END={end}"
    return {"chrom": chrom, "pos": pos, "id": rid, "ref": "N", "alt": alt,
            "info": info, "qual": str(qual), "filter": filt,
            "svtype": svtype, "svlen": svlen, "end": end}


# Delly: some PASS, one LowQual, one singleton, one DUP
RECORDS_COLLAPSE_DELLY = [
    col("chr1", 1000,  "DEL_D_01", "DEL", 1000, 2000,  qual=50,  filt="PASS"),
    col("chr1", 3000,  "DEL_D_02", "DEL", 1000, 4000,  qual=30,  filt="LowQual"),
    col("chr1", 5000,  "DEL_D_03", "DEL", 1000, 6000,  qual=40,  filt="PASS"),
    col("chr1", 9000,  "DUP_D_01", "DUP", 1000, 10000, qual=50,  filt="PASS"),
]

# Manta: all PASS, one singleton, same positions as Delly (some exact, some shifted)
RECORDS_COLLAPSE_MANTA = [
    col("chr1", 1000,  "DEL_M_01", "DEL", 1000, 2000,  qual=80,  filt="PASS"),
    col("chr1", 3100,  "DEL_M_02", "DEL", 1000, 4100,  qual=80,  filt="PASS"),
    col("chr1", 7000,  "DEL_M_03", "DEL", 1000, 8000,  qual=80,  filt="PASS"),
    col("chr1", 9000,  "DUP_M_01", "DUP", 1000, 10000, qual=80,  filt="PASS"),
]

# Multiallelic fixture: one record with two ALTs to test the error path.
MULTIALLELIC_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##contig=<ID=chr1,length=248956422>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"""

RECORDS_MULTIALLELIC = [
    {"chrom": "chr1", "pos": 1000, "id": "MULTI_01",
     "ref": "N", "alt": "<DEL>,<DUP>", "qual": ".", "filter": "PASS",
     "info": "SVTYPE=DEL;END=2000"},
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="tests/fixtures")
    args = parser.parse_args()

    out = args.output_dir
    os.makedirs(out, exist_ok=True)
    os.makedirs(os.path.join(out, "expected"), exist_ok=True)

    print("Writing fixtureA.vcf.gz ...", flush=True)
    write_vcf_gz(RECORDS_A, os.path.join(out, "fixtureA.vcf.gz"))

    print("Writing fixtureB.vcf.gz ...", flush=True)
    write_vcf_gz(RECORDS_B, os.path.join(out, "fixtureB.vcf.gz"))

    print("Writing fixtureA.bcf ...", flush=True)
    write_bcf(os.path.join(out, "fixtureA.vcf.gz"), os.path.join(out, "fixtureA.bcf"))

    print("Writing fixtureB.bcf ...", flush=True)
    write_bcf(os.path.join(out, "fixtureB.vcf.gz"), os.path.join(out, "fixtureB.bcf"))

    print("Writing fixtureDB.vcf.gz ...", flush=True)
    write_vcf_gz(RECORDS_DB, os.path.join(out, "fixtureDB.vcf.gz"), header=DB_HEADER)
    print("Writing fixtureDB.bcf ...", flush=True)
    write_bcf(os.path.join(out, "fixtureDB.vcf.gz"), os.path.join(out, "fixtureDB.bcf"))

    # Chromosome order from HEADER for stable sort.
    chrom_order = {c: i for i, c in enumerate(["chr1", "chr2", "chrX"])}
    def sort_key(row):
        # (chrom-order, svtype-string, posA, idA, posB)
        return (chrom_order.get(row[0], 999), row[7], int(row[1]), row[3], int(row[4]))

    # Expected TSVs for three threshold scenarios. Under the new CLI rule
    # exactly one of --min-overlap / --min-jaccard is active per run; the
    # SIMILARITY column reports that metric for interval rows. BND rows
    # always use slop=100 proximity and are merged into every scenario.
    scenarios = [
        ("default",      "overlap", 0.5),
        ("strict",       "overlap", 0.8),
        ("jaccard_only", "jaccard", 0.5),
    ]
    bnd_rows = compute_bnd_expected(RECORDS_A, RECORDS_B, slop=100)
    for name, metric, thr in scenarios:
        rows = compute_expected(RECORDS_A, RECORDS_B, metric, thr) + bnd_rows
        rows.sort(key=sort_key)
        tsv_path = os.path.join(out, "expected", f"expected_{name}.tsv")
        write_expected_tsv(rows, tsv_path)
        print(f"  {tsv_path}: {len(rows)} rows ({metric}>={thr}; +{len(bnd_rows)} BND)")

    # Self-mode expected output: --self --min-overlap 0.5 fixtureA.
    self_rows = (compute_self_expected(RECORDS_A, "overlap", 0.5) +
                 compute_self_bnd_expected(RECORDS_A, slop=100))
    self_rows.sort(key=sort_key)
    self_path = os.path.join(out, "expected", "expected_self.tsv")
    write_expected_tsv(self_rows, self_path)
    print(f"  {self_path}: {len(self_rows)} rows (--self, overlap>=0.5 + BND slop=100)")

    print("Writing collapse_delly.vcf.gz ...")
    write_vcf_gz(RECORDS_COLLAPSE_DELLY,
                 os.path.join(out, "collapse_delly.vcf.gz"),
                 header=COLLAPSE_HEADER)
    print("Writing collapse_manta.vcf.gz ...")
    write_vcf_gz(RECORDS_COLLAPSE_MANTA,
                 os.path.join(out, "collapse_manta.vcf.gz"),
                 header=COLLAPSE_HEADER)

    print("Writing collapse_multiallelic.vcf.gz ...")
    vcf_text = MULTIALLELIC_HEADER + "\n"
    for r in RECORDS_MULTIALLELIC:
        vcf_text += (f"{r['chrom']}\t{r['pos']}\t{r['id']}\t{r['ref']}\t{r['alt']}"
                     f"\t{r['qual']}\t{r['filter']}\t{r['info']}\n")
    multi_path = os.path.join(out, "collapse_multiallelic.vcf.gz")
    subprocess.run(
        ["bgzip", "-c"],
        input=vcf_text.encode(), check=True,
        stdout=open(multi_path, "wb"),
        capture_output=False,
    )
    subprocess.run(["bcftools", "index", multi_path], check=True, capture_output=True)

    print("Done.")


if __name__ == "__main__":
    main()
