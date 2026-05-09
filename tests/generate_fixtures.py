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

# Records are dicts so edge cases (missing fields, conflicting fields, missing IDs)
# can be expressed without contorting a fixed tuple shape. Each dict carries:
#   - VCF columns: chrom, pos, id, ref, alt, info (preformatted INFO string)
#   - Matching truth: svtype, svlen, end (None if record is intentionally
#     malformed or unmatchable; compute_expected then skips it)
def normal(chrom, pos, rid, ref, alt, svtype, svlen, end):
    return {
        "chrom": chrom, "pos": pos, "id": rid, "ref": ref, "alt": alt,
        "info": f"SVTYPE={svtype};SVLEN={svlen};END={end}",
        "svtype": svtype, "svlen": svlen, "end": end,
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
    # TC11 BND and INS — silently skipped by matcha (unsupported svtype)
    {"chrom": "chr1", "pos": 31000, "id": "BND_A_11", "ref": "N",
     "alt": "N[chr1:32000[", "info": "SVTYPE=BND",
     "svtype": None, "svlen": None, "end": None},
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
    # TC11 BND in B too (silently skipped)
    {"chrom": "chr1", "pos": 31000, "id": "BND_B_11", "ref": "N",
     "alt": "N[chr1:32000[", "info": "SVTYPE=BND",
     "svtype": None, "svlen": None, "end": None},
]


def recip_overlap(pA, eA, pB, eB):
    ovl = min(eA, eB) - max(pA, pB)
    if ovl <= 0:
        return 0.0
    return ovl / min(eA - pA, eB - pB)


def jaccard(pA, eA, pB, eB):
    ovl = min(eA, eB) - max(pA, pB)
    if ovl <= 0:
        return 0.0
    union = max(eA, eB) - min(pA, pB)
    return ovl / union


def records_to_vcf(records):
    lines = [HEADER]
    for r in records:
        lines.append(
            f"{r['chrom']}\t{r['pos']}\t{r['id']}\t{r['ref']}\t{r['alt']}"
            f"\t.\tPASS\t{r['info']}"
        )
    return "\n".join(lines) + "\n"


def write_vcf_gz(records, path):
    vcf_text = records_to_vcf(records)
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


def compute_expected(records_a, records_b, min_recip=0.0, min_jac=0.0):
    """Return list of output rows (as lists) for given thresholds.

    Records whose svtype/svlen/end are None (intentionally malformed or
    unsupported) are skipped — they will never appear in matcha's output.
    """
    def matchable(r):
        return r["svtype"] is not None and r["svlen"] is not None and r["end"] is not None

    # Build lookup of B records by (chrom, svtype)
    b_by_key = {}
    for rec in records_b:
        if not matchable(rec):
            continue
        key = (rec["chrom"], rec["svtype"])
        b_by_key.setdefault(key, []).append(rec)

    rows = []
    for rec_a in records_a:
        if not matchable(rec_a):
            continue
        chrom, posA, idA, svtype, endA = (
            rec_a["chrom"], rec_a["pos"], rec_a["id"], rec_a["svtype"], rec_a["end"],
        )
        key = (chrom, svtype)
        if key not in b_by_key:
            continue
        for rec_b in b_by_key[key]:
            posB, idB, endB = rec_b["pos"], rec_b["id"], rec_b["end"]
            ro = recip_overlap(posA, endA, posB, endB)
            jac = jaccard(posA, endA, posB, endB)
            if ro >= min_recip and jac >= min_jac:
                rows.append([
                    chrom, str(posA), str(endA), idA,
                    str(posB), str(endB), idB,
                    svtype,
                    f"{ro:.6f}", f"{jac:.6f}",
                ])
    return rows


def write_expected_tsv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")


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

    # Expected TSVs for three threshold scenarios
    scenarios = [
        ("default",      0.5,  0.0),   # --min-reciprocal-overlap 0.5
        ("strict",       0.8,  0.8),   # --min-reciprocal-overlap 0.8 --min-jaccard 0.8
        ("jaccard_only", 0.0,  0.5),   # --min-jaccard 0.5
    ]
    for name, min_recip, min_jac in scenarios:
        rows = compute_expected(RECORDS_A, RECORDS_B, min_recip, min_jac)
        tsv_path = os.path.join(out, "expected", f"expected_{name}.tsv")
        write_expected_tsv(rows, tsv_path)
        print(f"  {tsv_path}: {len(rows)} rows (min_recip={min_recip}, min_jac={min_jac})")

    print("Done.")


if __name__ == "__main__":
    main()
