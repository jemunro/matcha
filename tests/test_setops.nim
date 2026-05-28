## Integration tests for matcha intersect / setdiff.
## Requires the matcha binary (nimble build), bcftools, and test fixtures.
## Run from project root: nim c --hints:off -r tests/test_setops.nim
echo "--------------- Test SetOps ---------------"

import std/[algorithm, os, osproc, sequtils, sets, strutils]
import test_utils

const BinPath  = "./matcha"
const FixtureA = "tests/fixtures/fixtureA.vcf.gz"
const FixtureB = "tests/fixtures/fixtureB.vcf.gz"
const SampleFixture = "tests/fixtures/collapse_caller1_1sample.vcf.gz"

proc run(args: string): (string, int) =
  runMatcha(BinPath, args, getEnv("MATCHA_TEST_TIMEOUT", "30"))

proc runMerged(args: string): (string, int) =
  let t = getEnv("MATCHA_TEST_TIMEOUT", "30")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>&1")

proc dataLines(s: string): seq[string] =
  ## Non-header (non-'#') lines of a VCF body.
  for line in s.strip.splitLines:
    if line.len > 0 and line[0] != '#':
      result.add(line)

proc idsOf(vcfText: string): HashSet[string] =
  ## ID column (3rd, 0-based index 2) of each VCF data record.
  for line in dataLines(vcfText):
    let cols = line.split('\t')
    if cols.len >= 3: result.incl(cols[2])

proc bcftoolsIds(path: string): HashSet[string] =
  ## Ground-truth ID set of a fixture, read independently via bcftools.
  let (outp, code) = execCmdEx("bcftools view " & path & " 2>/dev/null")
  doAssert code == 0, "bcftools view failed for " & path
  result = idsOf(outp)

# T01 — binary available
timed("T01", "binary available"):
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp
      quit(1)
  doAssert fileExists(BinPath), "binary not found: " & BinPath

# T02 — intersect: default metric, valid VCF output
timed("T02", "intersect: exit 0, emits VCF with #CHROM header"):
  let (outp, code) = run("intersect " & FixtureA & " " & FixtureB)
  doAssert code == 0, "intersect should exit 0, got: " & outp
  doAssert "#CHROM" in outp, "output should contain a #CHROM header line"

# T03 — intersect IDs == the A-side IDs that `match` reports as matched
timed("T03", "intersect A B == matched A IDs from `match`"):
  let (mOut, mCode) = run("match --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert mCode == 0
  var matchedA = initHashSet[string]()
  for line in dataLines(mOut):
    let cols = line.split('\t')
    if cols.len >= 3: matchedA.incl(cols[2])  # ID_A column
  let (iOut, iCode) = run("intersect --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert iCode == 0, "exit " & $iCode & ": " & iOut
  let intersectIds = idsOf(iOut)
  doAssert intersectIds == matchedA,
    "intersect IDs differ from matched A IDs\nintersect: " &
    intersectIds.toSeq.sorted.join(",") & "\nmatched:   " & matchedA.toSeq.sorted.join(",")

# T04 — setdiff is the complement of intersect over all A records
timed("T04", "intersect ⊎ setdiff partitions A exactly"):
  let allA = bcftoolsIds(FixtureA)
  let (iOut, ic) = run("intersect --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  let (sOut, sc) = run("setdiff --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert ic == 0 and sc == 0, "exits: " & $ic & ", " & $sc
  let iIds = idsOf(iOut)
  let sIds = idsOf(sOut)
  doAssert (iIds * sIds).len == 0, "intersect and setdiff overlap: " &
    (iIds * sIds).toSeq.sorted.join(",")
  doAssert (iIds + sIds) == allA,
    "intersect ∪ setdiff != all A records (|A|=" & $allA.len &
    ", union=" & $((iIds + sIds).len) & ")"

# T05 — records are emitted verbatim: sample/GT columns preserved.
# intersect(F, F): every record self-matches, so all are kept unchanged.
timed("T05", "verbatim output preserves FORMAT/genotype columns"):
  let (outp, code) = run("intersect --min-overlap 0.5 " & SampleFixture & " " & SampleFixture)
  doAssert code == 0, "exit " & $code & ": " & outp
  let src = bcftoolsIds(SampleFixture)
  doAssert idsOf(outp) == src, "intersect(F,F) should keep every record of F"
  # Each data line must retain FORMAT (col 9) + at least one sample (col 10).
  for line in dataLines(outp):
    let cols = line.split('\t')
    doAssert cols.len >= 10, "record lost FORMAT/sample columns: " & line
    doAssert cols[8] == "GT", "FORMAT column not preserved: " & line
    doAssert cols[9] in ["0/0", "0/1", "1/1", "1/0", "./."],
      "genotype value not preserved verbatim: " & cols[9]

# T06 — setdiff(F, F) is empty (every record self-matches in B)
timed("T06", "setdiff(F,F) is empty"):
  let (outp, code) = run("setdiff --min-overlap 0.5 " & SampleFixture & " " & SampleFixture)
  doAssert code == 0, "exit " & $code & ": " & outp
  doAssert dataLines(outp).len == 0, "setdiff(F,F) should emit no records"

# T07 — stricter threshold keeps fewer-or-equal intersect records
timed("T07", "stricter --min-overlap shrinks the intersect set"):
  let (loose, lc) = run("intersect --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  let (tight, tc) = run("intersect --min-overlap 0.99 " & FixtureA & " " & FixtureB)
  doAssert lc == 0 and tc == 0
  let looseIds = idsOf(loose)
  let tightIds = idsOf(tight)
  doAssert tightIds <= looseIds, "tight intersect not a subset of loose intersect"
  doAssert tightIds.len <= looseIds.len

# T08 — BND/INS exercised: a matched BND is kept by intersect, an unmatched
# BND is kept by setdiff. (BND_A_11 matches BND_B_11; BND_A_22 has no mate.)
timed("T08", "intersect keeps matched BND; setdiff keeps unmatched BND"):
  let (iOut, ic) = run("intersect --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  let (sOut, sc) = run("setdiff --min-overlap 0.5 " & FixtureA & " " & FixtureB)
  doAssert ic == 0 and sc == 0
  doAssert "BND_A_11" in idsOf(iOut), "matched BND_A_11 should survive intersect"
  doAssert "BND_A_22" in idsOf(sOut), "unmatched BND_A_22 should survive setdiff"
  doAssert "INS_A_30_seq" in idsOf(iOut), "matched INS_A_30_seq should survive intersect"

# T09 — --threads 2 output equals --threads 1 output (sorted records)
timed("T09", "--threads 2 == --threads 1"):
  let (o1, c1) = run("intersect --min-overlap 0.5 --threads 1 " & FixtureA & " " & FixtureB)
  let (o2, c2) = run("intersect --min-overlap 0.5 --threads 2 " & FixtureA & " " & FixtureB)
  doAssert c1 == 0 and c2 == 0, "exits: " & $c1 & ", " & $c2
  doAssert dataLines(o1).sorted == dataLines(o2).sorted,
    "multi-thread intersect output differs from single-thread"

# T10 — -o writes file; --write-index produces a .csi for .vcf.gz
timed("T10", "-o + --write-index writes bgzipped output and CSI"):
  let tmpOut = getTempDir() / "matcha_setops_out.vcf.gz"
  defer:
    if fileExists(tmpOut): removeFile(tmpOut)
    if fileExists(tmpOut & ".csi"): removeFile(tmpOut & ".csi")
  let (outp, code) = run(
    "intersect --min-overlap 0.5 --write-index -o " & tmpOut & " " &
    FixtureA & " " & FixtureB)
  doAssert code == 0, "exit " & $code & ": " & outp
  doAssert fileExists(tmpOut), "output file not created"
  doAssert fileExists(tmpOut & ".csi"), "CSI index not created"

# T11 — error paths: wrong positional count and mutually exclusive thresholds
timed("T11", "wrong arg count / mutually exclusive metrics exit non-zero"):
  block:
    let (outp, code) = runMerged("intersect " & FixtureA)
    doAssert code != 0, "1 positional should fail"
    doAssert "2 input" in outp, "error should mention expected 2 inputs: " & outp
  block:
    let (outp, code) = runMerged(
      "setdiff --min-overlap 0.5 --min-jaccard 0.5 " & FixtureA & " " & FixtureB)
    doAssert code != 0, "both thresholds should fail"
    doAssert "mutually exclusive" in outp.toLowerAscii, "error should mention exclusion: " & outp
