## Tests for src/matcha/preproc.nim — VCF splitting and temp BCF writing.
## Run from project root: nim c --hints:off -r tests/test_preproc.nim
## Requires fixtures: run tests/generate_fixtures.py first.
echo "--------------- Test Preproc ---------------"

import std/[os, sets, strutils, tables, tempfiles]
import hts
import test_utils
import matcha/utils
import matcha/preproc
import matcha/log

# Silence the `[matcha preproc WARN] ...` lines that would otherwise interleave
# with PASS/FAIL output for every test that exercises preprocessing.
setQuiet(true)

const FixtureA = "tests/fixtures/fixtureA.vcf.gz"
const FixtureB = "tests/fixtures/fixtureB.vcf.gz"

# Convenience: paths key for "this svtype, bin 0" — the bin all the 1000bp
# fixture records land in.
template delBin0(pp: PreprocOutput): string = pp.paths[(svDEL, 0)]

# A MatchConfig with --min-overlap 0.5 (same as smoke tests). buildWorkQueue
# uses the threshold for adjacent-bin pruning.
proc baseCfg(): MatchConfig =
  MatchConfig(metric: mOverlap, threshold: 0.5, nThreads: 1)

# P01 — preprocessVcf produces the expected (svtype, bin) keys + chroms
timed("P01", "preprocessVcf: (svtype, bin) keys + populated"):
  doAssert fileExists(FixtureA), "fixture missing: " & FixtureA
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  doAssert (svDEL, 0) in pp.paths, "missing (DEL, bin 0)"
  doAssert (svDUP, 1) in pp.paths, "missing (DUP, bin 1)"
  doAssert (svINV, 1) in pp.paths, "missing (INV, bin 1)"
  doAssert "chr1" in pp.populated.getOrDefault((svDEL, 0)), "DEL/bin0 missing chr1"
  doAssert "chr2" in pp.populated.getOrDefault((svDUP, 1)), "DUP/bin1 missing chr2"
  doAssert "chrX" in pp.populated.getOrDefault((svINV, 1)), "INV/bin1 missing chrX"

# P02 — BND records are kept; INS and TRA are skipped.
timed("P02", "preprocessVcf: BND kept; INS and TRA excluded"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  doAssert (svBND, 0) in pp.paths,
    "BND records should be in the (svBND, bin 0) temp BCF"
  for key in pp.paths.keys:
    doAssert key.svtype != svINS, "INS should be excluded"
    doAssert key.svtype != svTRA, "TRA should be excluded"

# P02b — BND slim record carries authoritative CHR2/POS2 from ALT parse.
timed("P03", "preprocessVcf: BND slim record has CHR2/POS2 parsed from ALT"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  doAssert (svBND, 0) in pp.paths, "no BND temp BCF"
  let path = pp.paths[(svBND, 0)]
  var vcf: VCF
  doAssert open(vcf, path), "cannot open BND temp BCF: " & path
  var pos2Data: seq[int32]
  var chr2Data: string
  var foundA11 = false
  for v in vcf:
    if $v.ID == "BND_A_11":
      doAssert v.info().get("CHR2", chr2Data) == Status.OK, "CHR2 missing on BND_A_11"
      doAssert chr2Data == "chr1", "CHR2 should be chr1, got " & chr2Data
      doAssert v.info().get("POS2", pos2Data) == Status.OK, "POS2 missing on BND_A_11"
      doAssert pos2Data.len > 0 and int64(pos2Data[0]) == 32000,
        "POS2 should be 32000, got " & $pos2Data
      foundA11 = true
      break
  vcf.close()
  doAssert foundA11, "BND_A_11 missing from BND temp BCF"

# P03 — extracted fields are correct (POS, END, ID)
timed("P04", "preprocessVcf: DEL_A_01 has correct POS=1000, END=2000"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  let path = pp.delBin0
  var vcf: VCF
  doAssert open(vcf, path), "cannot open temp BCF: " & path
  var found = false
  var endData: seq[int32]
  for v in vcf:
    if $v.ID == "DEL_A_01":
      doAssert v.POS == 1000, "POS mismatch for DEL_A_01: " & $v.POS
      if v.info().get("END", endData) == Status.OK and endData.len > 0:
        doAssert endData[0] == 2000, "END mismatch for DEL_A_01: " & $endData[0]
      found = true
      break
  vcf.close()
  doAssert found, "DEL_A_01 not found in temp BCF"

# P04 — every temp BCF has a CSI index
timed("P05", "preprocessVcf: every temp BCF has a .csi index file"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  for key, path in pp.paths:
    doAssert fileExists(path & ".csi"),
      "missing CSI index for " & $key.svtype & "/bin" & $key.bin & ": " & path

# P05 — every emitted job key is in both A and B
timed("P06", "buildWorkQueue: every job (chrom, svtype, binA) is reachable in both"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let ppA = preprocessVcf(FixtureA, tmpDir, "A")
  let ppB = preprocessVcf(FixtureB, tmpDir, "B")
  let (jobs, _) = buildWorkQueue(ppA, ppB, baseCfg())
  doAssert jobs.len > 0, "expected at least one job"
  for job in jobs:
    doAssert (job.svtype, job.binA) in ppA.paths,
      "job's (svtype, binA) not in A: " & $job.svtype & "/bin" & $job.binA
    doAssert job.chrom in ppA.populated.getOrDefault((job.svtype, job.binA)),
      "chrom not in A's populated set: " & job.chrom & "/" & $job.svtype & "/bin" & $job.binA
    block chromInB:
      for binB in job.binsB.keys:
        if job.chrom in ppB.populated.getOrDefault((job.svtype, binB)):
          break chromInB
      doAssert false, "chrom not in any adjacent B bin: " & job.chrom & "/" & $job.svtype
    doAssert job.binsB.len > 0, "job has no adjacent B bins"
    for binB in job.binsB.keys:
      doAssert (job.svtype, binB) in ppB.paths,
        "job binB references missing B path"

# P06 — chrom order in jobs respects A's input header order
timed("P07", "buildWorkQueue: jobs ordered by VCF header chrom order"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let ppA = preprocessVcf(FixtureA, tmpDir, "A")
  let ppB = preprocessVcf(FixtureB, tmpDir, "B")
  let (jobs, _) = buildWorkQueue(ppA, ppB, baseCfg())
  # Map chrom → index in A's header order; jobs' chroms must be monotonic.
  var idx: Table[string, int]
  for i, c in ppA.chromOrder: idx[c] = i
  var prev = -1
  for job in jobs:
    let cur = idx.getOrDefault(job.chrom, high(int))
    doAssert cur >= prev,
      "chrom order regressed at job " & job.chrom &
      " (prev=" & $prev & " cur=" & $cur & ")"
    prev = cur

# P07 — preprocessVcf works on .bcf inputs
timed("P08", "preprocessVcf: accepts .bcf input"):
  const FixtureA_bcf = "tests/fixtures/fixtureA.bcf"
  doAssert fileExists(FixtureA_bcf), "fixture missing: " & FixtureA_bcf
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA_bcf, tmpDir, "A_bcf")
  doAssert (svDEL, 0) in pp.paths, "missing (DEL, bin 0) from .bcf input"
  doAssert "chr1" in pp.populated.getOrDefault((svDEL, 0)), "DEL/bin0/chr1 should be present"

# Helper: collect every record in a per-(svtype, bin) BCF.
proc readRecords(path: string): seq[tuple[id: string, pos: int64, endPos: int64,
                                          svtype: string, infoNames: seq[string],
                                          svlen: int64]] =
  var vcf: VCF
  doAssert open(vcf, path), "cannot open " & path
  var endData, svlenData: seq[int32]
  var svtypeStr: string
  for v in vcf:
    var infoNames: seq[string]
    for fld in v.info.fields:
      infoNames.add(fld.name)
    var endVal: int64 = 0
    if v.info.get("END", endData) == Status.OK and endData.len > 0:
      endVal = int64(endData[0])
    var svlenVal: int64 = 0
    if v.info.get("SVLEN", svlenData) == Status.OK and svlenData.len > 0:
      svlenVal = int64(svlenData[0])
    var sv: string
    if v.info.get("SVTYPE", svtypeStr) == Status.OK and svtypeStr.len > 0:
      sv = svtypeStr
    result.add((id: $v.ID, pos: v.POS, endPos: endVal, svtype: sv,
                infoNames: infoNames, svlen: svlenVal))
  vcf.close()

# P08 — INFO is slimmed to per-SVTYPE keep-set only after preprocessing
timed("P09", "preprocessVcf: temp BCF records have only keep-set INFO fields"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  let allowedInterval = ["END", "SRC_INDEX"].toHashSet
  let allowedBnd      = ["CHR2", "POS2", "SRC_INDEX"].toHashSet
  for key, path in pp.paths:
    let allowed = if key.svtype == svBND: allowedBnd else: allowedInterval
    for rec in readRecords(path):
      for name in rec.infoNames:
        doAssert name in allowed,
          "unexpected INFO field " & name & " in " & rec.id &
          " (svtype=" & $key.svtype & ")"

# P09 — record with ID="." gets synthesized (CHROM_POS_SVTYPE_LINENUMBER)
timed("P10", "preprocessVcf: missing ID is synthesized"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  var found = false
  for rec in readRecords(pp.delBin0):
    # TC14 record sits at chr1:37000 in the input; synthetic ID format is
    # CHROM_POS_SVTYPE_LINENUMBER. The lineno is the 1-based input order.
    if rec.id.startsWith("chr1_37000_DEL_"):
      found = true
      break
  doAssert found, "expected synthesized ID prefix chr1_37000_DEL_<lineno>"

# P10 — symbolic ALT only (no INFO/SVTYPE) is routed to the right SVTYPE BCF
timed("P11", "preprocessVcf: ALT-symbolic SVTYPE (TC12) routed to (DEL, 0)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  var found = false
  for rec in readRecords(pp.delBin0):
    if rec.id == "DEL_A_12_alt_only":
      found = true
      break
  doAssert found, "DEL_A_12_alt_only missing from (DEL, 0) BCF"

# P11 — SVTYPE conflict (INFO=DUP, ALT=<DEL>) → ALT wins; record lands in DEL BCF
timed("P12", "preprocessVcf: ALT wins on SVTYPE conflict (TC13)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  var inDel = false
  for rec in readRecords(pp.delBin0):
    if rec.id == "DEL_A_13_conflict":
      inDel = true
      break
  doAssert inDel, "DEL_A_13_conflict should be in (DEL, 0) BCF"
  for key, path in pp.paths:
    if key.svtype == svDUP:
      for rec in readRecords(path):
        doAssert rec.id != "DEL_A_13_conflict", "should not be in any DUP bin"

# P12 — record with neither END nor SVLEN is skipped
timed("P13", "preprocessVcf: missing END and SVLEN → skipped (TC15)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  for path in pp.paths.values:
    for rec in readRecords(path):
      doAssert rec.id != "DEL_A_15_no_end_no_svlen",
        "TC15 should have been skipped (no END, no SVLEN)"

# P13 — record with END < POS is skipped
timed("P14", "preprocessVcf: END <= POS → skipped (TC16)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  for path in pp.paths.values:
    for rec in readRecords(path):
      doAssert rec.id != "DEL_A_16_bad_end",
        "TC16 should have been skipped (END < POS)"

# P15 — SRC_INDEX is present on every slim BCF record and values are unique
timed("P15", "preprocessVcf: SRC_INDEX present and unique on slim BCF records"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")

  var slim: VCF
  doAssert open(slim, pp.delBin0), "cannot open slim DEL BCF"
  var idxData: seq[int32]
  var seen: seq[int32]
  for v in slim:
    doAssert v.info.get("SRC_INDEX", idxData) == Status.OK,
      "SRC_INDEX missing on record " & $v.ID
    doAssert idxData.len == 1,
      "SRC_INDEX should be Number=1, got len=" & $idxData.len
    doAssert idxData[0] >= 0,
      "SRC_INDEX should be non-negative, got " & $idxData[0]
    doAssert idxData[0] notin seen,
      "SRC_INDEX " & $idxData[0] & " duplicated"
    seen.add(idxData[0])
  slim.close()
  doAssert seen.len > 0, "no records in slim DEL BCF"

# P16 — large SV lands in a non-zero bin (DEL_A_08 = 5000bp → bin 3)
timed("P16", "preprocessVcf: 5000bp DEL_A_08 lands in (DEL, bin 3)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  doAssert (svDEL, 3) in pp.paths,
    "expected (DEL, bin 3) populated for DEL_A_08"
  doAssert (svDEL, 3) in pp.populated,
    "bin 3 should be in populated for DEL"
  var foundLarge = false
  for rec in readRecords(pp.paths[(svDEL, 3)]):
    if rec.id == "DEL_A_08":
      doAssert rec.pos == 17000, "DEL_A_08 POS"
      doAssert rec.endPos == 22000, "DEL_A_08 END"
      foundLarge = true
      break
  doAssert foundLarge, "DEL_A_08 missing from (DEL, 3) BCF"
  # And it should NOT be in (DEL, 0).
  for rec in readRecords(pp.delBin0):
    doAssert rec.id != "DEL_A_08", "DEL_A_08 should not be in (DEL, 0)"
