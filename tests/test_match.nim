## Integration tests for matcha match.
## Requires the matcha binary (nimble build) and test fixtures (generate_fixtures.py).
## Run from project root: nim c --hints:off -r tests/test_match.nim
echo "--------------- Test Match ---------------"

import std/[algorithm, os, osproc, sequtils, sets, strutils]
import test_utils

const BinPath  = "./matcha"
const FixtureA = "tests/fixtures/fixtureA.vcf.gz"
const FixtureB = "tests/fixtures/fixtureB.vcf.gz"

proc run(args: string): (string, int) =
  ## Run matcha and capture stdout only. Preproc warnings on stderr would
  ## otherwise pollute parseTsv and break column-shape assertions.
  let t = getEnv("MATCHA_TEST_TIMEOUT", "30")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>/dev/null")

proc runMerged(args: string): (string, int) =
  ## Run matcha and capture stdout+stderr (for tests that assert on errors).
  let t = getEnv("MATCHA_TEST_TIMEOUT", "30")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>&1")

proc parseTsv(s: string): seq[seq[string]] =
  ## Skip blank lines and any line starting with '#' (header / comments).
  for line in s.strip.splitLines:
    if line.len > 0 and line[0] != '#':
      result.add(line.split('\t'))

# M01 — binary available
timed("M01", "binary available"):
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp
      quit(1)
  doAssert fileExists(BinPath), "binary not found: " & BinPath

# M02 — no threshold exits non-zero with informative error (on stderr)
timed("M02", "no threshold flag: exits non-zero, mentions flag names"):
  let (outp, code) = runMerged("match " & FixtureA & " " & FixtureB)
  doAssert code != 0, "should exit non-zero without threshold"
  doAssert "min-overlap" in outp or "min-jaccard" in outp,
    "error should mention threshold flags, got: " & outp

# M02b — passing BOTH --min-overlap and --min-jaccard is an error (xor).
timed("M02b", "both thresholds: mutually exclusive error"):
  let (outp, code) = runMerged(
    "match --min-overlap 0.5 --min-jaccard 0.5 " & FixtureA & " " & FixtureB)
  doAssert code != 0, "passing both should fail"
  doAssert "mutually exclusive" in outp.toLowerAscii or
           "exactly one" in outp.toLowerAscii,
    "error should mention mutual exclusion, got: " & outp

# M03 — TC01 exact match: SIMILARITY=1.0
timed("M03", "exact match DEL_A_01/DEL_B_01: similarity = 1.0"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0, "exit " & $code & ": " & outp
  let rows = parseTsv(outp)
  var found = false
  for row in rows:
    if row.len >= 9 and row[3] == "DEL_A_01" and row[6] == "DEL_B_01":
      doAssert abs(parseFloat(row[8]) - 1.0) < 1e-5, "SIMILARITY != 1.0: " & row[8]
      found = true
  doAssert found, "DEL_A_01/DEL_B_01 pair not found in output"

# M04 — TC02 partial overlap above threshold is emitted
timed("M04", "partial overlap above threshold: DEL_A_02 emitted"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  var found = false
  for row in rows:
    if row.len >= 9 and row[3] == "DEL_A_02":
      found = true
  doAssert found, "DEL_A_02 partial overlap should be emitted at threshold 0.5"

# M05 — TC03 below threshold NOT emitted
timed("M05", "below-threshold overlap: DEL_A_03/DEL_B_03 not emitted"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  for row in rows:
    if row.len >= 9:
      doAssert not (row[3] == "DEL_A_03" and row[6] == "DEL_B_03"),
        "DEL_A_03/DEL_B_03 overlap=0.4 should be filtered at threshold 0.5"

# M06 — TC05 SVTYPE mismatch produces no output
timed("M06", "SVTYPE mismatch DEL_A_05/DUP_B_05: no match emitted"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  for row in rows:
    if row.len >= 9:
      doAssert not (row[3] == "DEL_A_05" and row[6] == "DUP_B_05"),
        "DEL/DUP SVTYPE mismatch should not match"

# M07 — TC06 multiple matches: both B records emitted for DEL_A_06
timed("M07", "multiple matches: 2 B records match DEL_A_06"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  var matchCount = 0
  for row in rows:
    if row.len >= 9 and row[3] == "DEL_A_06":
      inc matchCount
  doAssert matchCount == 2, "expected 2 matches for DEL_A_06, got " & $matchCount

# M08 — TC07 unmatched A record: DEL_A_07 absent from output
timed("M08", "unmatched A record: DEL_A_07 not in output"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  for row in rows:
    if row.len >= 9:
      doAssert row[3] != "DEL_A_07", "DEL_A_07 should have no match in B"

# M09 — TC09 multi-chromosome: chr2 and chrX rows present
timed("M09", "multi-chromosome: chr2 and chrX results present"):
  let (outp, code) = run("match --min-overlap 0.9 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  var hasChr2 = false
  var hasChrX = false
  for row in rows:
    if row.len >= 9:
      if row[0] == "chr2": hasChr2 = true
      if row[0] == "chrX": hasChrX = true
  doAssert hasChr2, "chr2 results missing"
  doAssert hasChrX, "chrX results missing"

# M10 — --min-jaccard alone: TC08 DEL_A_08 excluded (jaccard=0.2 < 0.5)
timed("M10", "--min-jaccard 0.5 alone: DEL_A_08 excluded (jaccard=0.2)"):
  let (outp, code) = run("match --min-jaccard 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  for row in rows:
    if row.len >= 9:
      doAssert row[3] != "DEL_A_08",
        "DEL_A_08 has jaccard=0.2 and should be excluded at --min-jaccard 0.5"

# M11 — output has 9 tab-separated columns, similarity in [0,1].
# Interval rows have numeric END_A/END_B; BND rows emit ".".
timed("M11", "output: 9 columns, similarity in [0,1]"):
  let (outp, code) = run("match --min-overlap 0.9 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  let rows = parseTsv(outp)
  doAssert rows.len > 0, "expected output rows"
  for row in rows:
    doAssert row.len == 9, "expected 9 columns, got " & $row.len & ": " & row.join("\t")
    doAssert row[1].parseInt > 0, "POS_A must be positive integer"
    doAssert row[4].parseInt > 0, "POS_B must be positive integer"
    let sim = parseFloat(row[8])
    doAssert sim >= 0.0 and sim <= 1.0, "SIMILARITY out of [0,1]: " & $sim
    if row[7] == "BND":
      doAssert row[2] == "." and row[5] == ".", "BND END cols should be '.'"
    else:
      doAssert row[2].parseInt > 0 and row[5].parseInt > 0,
        "interval END cols should be positive"

# M12 — --output writes to file
timed("M12", "--output: results written to file"):
  let tmpOut = getTempDir() / "matcha_test_out.tsv"
  defer: (if fileExists(tmpOut): removeFile(tmpOut))
  let (outp, code) = run(
    "match --min-overlap 0.9 --output " & tmpOut & " " &
    FixtureA & " " & FixtureB)
  doAssert code == 0, "exit " & $code & ": " & outp
  doAssert fileExists(tmpOut), "output file not created"
  doAssert readFile(tmpOut).len > 0, "output file is empty"

# M13 — --threads 2 output (sorted) equals --threads 1 output (sorted)
timed("M13", "--threads 2 produces same output as --threads 1"):
  let (outp1, c1) = run("match --min-overlap 0.5 --threads 1 " & FixtureA & " " & FixtureB)
  let (outp2, c2) = run("match --min-overlap 0.5 --threads 2 " & FixtureA & " " & FixtureB)
  doAssert c1 == 0 and c2 == 0,
    "exits: " & $c1 & ", " & $c2 & "\n" & outp1 & "\n" & outp2
  let sorted1 = outp1.strip.splitLines.sorted.join("\n")
  let sorted2 = outp2.strip.splitLines.sorted.join("\n")
  doAssert sorted1 == sorted2, "multi-thread output differs from single-thread"

# M14 — output begins with ##matcha_metric= line, then a #-prefixed header
# with the SIMILARITY column.
timed("M14", "output: ##matcha_metric preamble + SIMILARITY header"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0, "exit " & $code & ": " & outp
  let lines = outp.strip.splitLines
  doAssert lines.len > 1, "expected preamble + header + rows"
  doAssert lines[0] == "##matcha_metric=overlap",
    "first line should be ##matcha_metric=overlap, got: " & lines[0]
  let header = lines[1]
  doAssert header.startsWith("#CHROM"), "second line should be the #CHROM header"
  doAssert "SIMILARITY" in header, "header missing SIMILARITY column: " & header
  doAssert "OVERLAP" notin header, "OVERLAP column should be gone"
  doAssert "JACCARD" notin header, "JACCARD column should be gone"

# M14b — --min-jaccard switches the metric line accordingly.
timed("M14b", "##matcha_metric=jaccard when --min-jaccard is active"):
  let (outp, code) = run("match --min-jaccard 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0, "exit " & $code & ": " & outp
  let lines = outp.strip.splitLines
  doAssert lines[0] == "##matcha_metric=jaccard",
    "got: " & lines[0]

# M15 — .bcf inputs produce identical output to .vcf.gz inputs
timed("M15", ".bcf inputs: same output as .vcf.gz inputs"):
  const FixtureA_bcf = "tests/fixtures/fixtureA.bcf"
  const FixtureB_bcf = "tests/fixtures/fixtureB.bcf"
  doAssert fileExists(FixtureA_bcf), "fixture missing: " & FixtureA_bcf
  doAssert fileExists(FixtureB_bcf), "fixture missing: " & FixtureB_bcf
  let (vOut, vCode) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  let (bOut, bCode) = run("match --min-overlap 0.5 " & FixtureA_bcf & " " & FixtureB_bcf)
  doAssert vCode == 0 and bCode == 0,
    "exits: vcf=" & $vCode & " bcf=" & $bCode & "\n" & vOut & "\n" & bOut
  let vSorted = vOut.strip.splitLines.sorted.join("\n")
  let bSorted = bOut.strip.splitLines.sorted.join("\n")
  doAssert vSorted == bSorted, "BCF output differs from VCF output"

# ---------------------------------------------------------------------------
# Self-mode tests (S-prefix)
# ---------------------------------------------------------------------------

const ExpectedSelf = "tests/fixtures/expected/expected_self.tsv"

# S01 — --self emits exactly the expected pair set against fixtureA.
timed("S01", "--self matches expected_self.tsv pair set"):
  doAssert fileExists(ExpectedSelf), "fixture missing: " & ExpectedSelf
  let (outp, code) = run("match --self --min-overlap 0.5 " & FixtureA)
  doAssert code == 0, "exit " & $code & ": " & outp
  let actual = outp.strip.splitLines.filterIt(it.len > 0 and it[0] != '#').sorted.join("\n")
  let expected = readFile(ExpectedSelf).strip.splitLines.sorted.join("\n")
  doAssert actual == expected,
    "self-mode output differs from expected\n--- actual ---\n" & actual &
    "\n--- expected ---\n" & expected

# S02 — --self never emits a self-self pair (ID_A == ID_B)
timed("S02", "--self: no row has ID_A == ID_B"):
  let (outp, code) = run("match --self --min-overlap 0.5 " & FixtureA)
  doAssert code == 0
  for row in parseTsv(outp):
    if row.len >= 10:
      doAssert row[3] != row[6], "self-self pair leaked: " & row.join("\t")

# S03 — --self emits each unordered (X,Y) pair at most once
timed("S03", "--self: no duplicate symmetric pairs"):
  let (outp, code) = run("match --self --min-overlap 0.5 " & FixtureA)
  doAssert code == 0
  var seen: HashSet[string]
  for row in parseTsv(outp):
    if row.len >= 10:
      # canonical pair key: lexicographically min/max of the two IDs
      let key = (if row[3] < row[6]: row[3] & "|" & row[6]
                 else: row[6] & "|" & row[3])
      doAssert key notin seen, "duplicate pair: " & key
      seen.incl(key)

# S04 — --self with one positional input succeeds
timed("S04", "--self: 1 positional arg is accepted"):
  let (outp, code) = run("match --self --min-overlap 0.5 " & FixtureA)
  doAssert code == 0, "expected exit 0, got " & $code & ": " & outp

# S05 — --self with 0 or 2 positionals is an error mentioning --self
timed("S05", "--self: 0 or 2 positionals errors out, mentions --self"):
  block:
    let (outp, code) = runMerged("match --self --min-overlap 0.5")
    doAssert code != 0, "expected non-zero exit with 0 positionals"
    doAssert "--self" in outp, "error message should mention --self: " & outp
  block:
    let (outp, code) = runMerged(
      "match --self --min-overlap 0.5 " & FixtureA & " " & FixtureB)
    doAssert code != 0, "expected non-zero exit with 2 positionals in --self"
    doAssert "--self" in outp, "error message should mention --self: " & outp

# S06 — --self --threads 2 is a permutation of --threads 1 output
timed("S06", "--self --threads 2 produces same output as --threads 1"):
  let (out1, c1) = run("match --self --min-overlap 0.5 --threads 1 " & FixtureA)
  let (out2, c2) = run("match --self --min-overlap 0.5 --threads 2 " & FixtureA)
  doAssert c1 == 0 and c2 == 0, "exits: " & $c1 & ", " & $c2
  let s1 = out1.strip.splitLines.sorted.join("\n")
  let s2 = out2.strip.splitLines.sorted.join("\n")
  doAssert s1 == s2, "multi-thread self output differs from single-thread"

# S07 — .bcf input parity in --self mode
timed("S07", "--self: .bcf input matches .vcf.gz input"):
  const FixtureA_bcf = "tests/fixtures/fixtureA.bcf"
  let (vOut, vc) = run("match --self --min-overlap 0.5 " & FixtureA)
  let (bOut, bc) = run("match --self --min-overlap 0.5 " & FixtureA_bcf)
  doAssert vc == 0 and bc == 0
  let vs = vOut.strip.splitLines.sorted.join("\n")
  let bs = bOut.strip.splitLines.sorted.join("\n")
  doAssert vs == bs, "BCF self-mode output differs from VCF"

# ---------------------------------------------------------------------------
# BND tests (B-prefix)
# ---------------------------------------------------------------------------

# B01 — BND_A_11 / BND_B_11 identical intra-chrom mate → sim=1.0.
timed("B01", "BND identical mate pair: similarity=1.0"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0, "exit " & $code & ": " & outp
  var found = false
  for row in parseTsv(outp):
    if row.len >= 9 and row[3] == "BND_A_11" and row[6] == "BND_B_11":
      doAssert row[7] == "BND", "expected SVTYPE=BND, got " & row[7]
      doAssert row[2] == "." and row[5] == ".", "BND END_A/END_B must be '.'"
      doAssert abs(parseFloat(row[8]) - 1.0) < 1e-5,
        "BND_A_11/B_11 sim != 1.0: " & row[8]
      found = true
  doAssert found, "BND_A_11/BND_B_11 missing from output"

# B02 — BND_A_20 / BND_B_20 inter-chrom identical: sim=1.0.
timed("B02", "BND inter-chrom identical: sim=1.0"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  var found = false
  for row in parseTsv(outp):
    if row.len >= 9 and row[3] == "BND_A_20" and row[6] == "BND_B_20":
      doAssert abs(parseFloat(row[8]) - 1.0) < 1e-5
      found = true
  doAssert found, "BND_A_20/BND_B_20 missing"

# B03 — BND_A_21 / BND_B_21 offsets: sim=0.6 default slop.
timed("B03", "BND offset 50/30: sim=0.60 at default slop 100"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  var found = false
  for row in parseTsv(outp):
    if row.len >= 9 and row[3] == "BND_A_21" and row[6] == "BND_B_21":
      doAssert abs(parseFloat(row[8]) - 0.6) < 1e-5,
        "BND_A_21/B_21 sim != 0.6: " & row[8]
      found = true
  doAssert found, "BND_A_21/BND_B_21 missing"

# B04 — Strict slop rejects BND_A_21 / BND_B_21 (dPOS=50 >= slop=20).
timed("B04", "--bnd-slop 20 rejects offset-50 pair"):
  let (outp, code) = run("match --min-overlap 0.5 --bnd-slop 20 " &
                         FixtureA & " " & FixtureB)
  doAssert code == 0
  for row in parseTsv(outp):
    if row.len >= 9:
      doAssert not (row[3] == "BND_A_21" and row[6] == "BND_B_21"),
        "BND_A_21/B_21 should be rejected at slop=20"
  # But BND_A_20/BND_B_20 (offsets 0/0) should still match.
  var foundExact = false
  for row in parseTsv(outp):
    if row.len >= 9 and row[3] == "BND_A_20" and row[6] == "BND_B_20":
      foundExact = true
  doAssert foundExact, "exact BND match should survive at any positive slop"

# B05 — CHR2 mismatch (BND_A_23 chrX vs BND_B_23 chr2): no row.
timed("B05", "CHR2 mismatch: BND_A_23 / BND_B_23 not paired"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  for row in parseTsv(outp):
    if row.len >= 9:
      doAssert not (row[3] == "BND_A_23" and row[6] == "BND_B_23"),
        "CHR2-mismatched pair should not be emitted"

# B06 — BND_A_22 has no B mate at chr1:80000 → absent.
timed("B06", "unmatched BND: BND_A_22 not in output"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  for row in parseTsv(outp):
    if row.len >= 9:
      doAssert row[3] != "BND_A_22", "BND_A_22 should be unmatched"

# B07 — Malformed BND ALT and TRA records are dropped at preproc
#       (warn-skip on stderr); the merged stderr should mention them.
timed("B07", "preproc warn-skips malformed BND and TRA on stderr"):
  let (outp, _) = runMerged(
    "match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert "malformed_bnd" in outp or "BND_A_24" in outp,
    "expected warn-skip for malformed BND_A_24, got stderr:\n" & outp
  doAssert "unsupported_tra" in outp or "TRA_A_25" in outp,
    "expected warn-skip for TRA_A_25, got stderr:\n" & outp

# B08 — --self pairs BND_A_27a / BND_A_27b once (sim=0.80, aOff<bOff).
timed("B08", "--self: BND mate pair appears exactly once"):
  let (outp, code) = run("match --self --min-overlap 0.5 " & FixtureA)
  doAssert code == 0, "exit " & $code & ": " & outp
  var count = 0
  var simVal = 0.0
  for row in parseTsv(outp):
    if row.len >= 9 and
       ((row[3] == "BND_A_27a" and row[6] == "BND_A_27b") or
        (row[3] == "BND_A_27b" and row[6] == "BND_A_27a")):
      inc count
      simVal = parseFloat(row[8])
  doAssert count == 1, "expected 1 row for BND_A_27a/27b, got " & $count
  doAssert abs(simVal - 0.8) < 1e-5, "self-pair sim != 0.8: " & $simVal

# B09 — Mixed DEL+BND callset: both DEL and BND rows appear in one run.
timed("B09", "mixed DEL+BND run produces both kinds of rows"):
  let (outp, code) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert code == 0
  var sawDel = false
  var sawBnd = false
  for row in parseTsv(outp):
    if row.len >= 9:
      if row[7] == "DEL": sawDel = true
      if row[7] == "BND": sawBnd = true
  doAssert sawDel and sawBnd,
    "expected both DEL and BND rows; saw DEL=" & $sawDel & " BND=" & $sawBnd
