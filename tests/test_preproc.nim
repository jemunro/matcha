## Tests for src/matcha/preproc.nim — VCF splitting and temp BCF writing.
## Run from project root: nim c --hints:off -r tests/test_preproc.nim
## Requires fixtures: run tests/generate_fixtures.py first.
echo "--------------- Test Preproc ---------------"

import std/[os, sets, strutils, tables, tempfiles]
import hts
import hts/private/hts_concat
import test_utils
import matcha/utils
import matcha/preproc

const FixtureA = "tests/fixtures/fixtureA.vcf.gz"
const FixtureB = "tests/fixtures/fixtureB.vcf.gz"

# Local bindings for the P15 round-trip test. bgzf_seek is not bound by
# vendored hts-nim; the parallel struct mirrors hts-nim's VCF prefix so we
# can reach its private htsFile* / BGZF* without patching the vendored copy.
proc bgzf_seek(fp: ptr BGZF, pos: int64, whence: cint): int64
  {.cdecl, importc: "bgzf_seek", dynlib: "libhts.so".}

type VcfPrivT = ref object of RootObj
  hts: ptr htsFile

proc rawHts(v: VCF): ptr htsFile {.inline.} = cast[VcfPrivT](v).hts
proc rawBgzf(v: VCF): ptr BGZF {.inline.} = cast[VcfPrivT](v).hts.fp.bgzf

# P01 — preprocessVcf produces the expected SVTYPE keys with the right chroms
timed("P01", "preprocessVcf: SVTYPE keys + chromsBySvtype populated for A"):
  doAssert fileExists(FixtureA), "fixture missing: " & FixtureA
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  doAssert svDEL in pp.paths, "missing DEL"
  doAssert svDUP in pp.paths, "missing DUP"
  doAssert svINV in pp.paths, "missing INV"
  doAssert "chr1" in pp.chromsBySvtype[svDEL], "DEL should include chr1"
  doAssert "chr1" in pp.chromsBySvtype[svDUP], "DUP should include chr1"
  doAssert "chr2" in pp.chromsBySvtype[svDUP], "DUP should include chr2"
  doAssert "chr1" in pp.chromsBySvtype[svINV], "INV should include chr1"
  doAssert "chrX" in pp.chromsBySvtype[svINV], "INV should include chrX"

# P02 — BND and INS records are excluded
timed("P02", "preprocessVcf: BND and INS excluded from paths"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  doAssert svBND notin pp.paths, "BND should be excluded"
  doAssert svINS notin pp.paths, "INS should be excluded"

# P03 — extracted fields are correct (POS, END, ID)
timed("P03", "preprocessVcf: DEL_A_01 has correct POS=1000, END=2000"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  let path = pp.paths[svDEL]
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
timed("P04", "preprocessVcf: every temp BCF has a .csi index file"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  for svt, path in pp.paths:
    doAssert fileExists(path & ".csi"),
      "missing CSI index for " & $svt & ": " & path

# P05 — work queue keys are a subset of both A and B
timed("P05", "buildWorkQueue: every job (chrom,svtype) is in both A and B"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let ppA = preprocessVcf(FixtureA, tmpDir, "A")
  let ppB = preprocessVcf(FixtureB, tmpDir, "B")
  let jobs = buildWorkQueue(ppA, ppB, tmpDir)
  doAssert jobs.len > 0, "expected at least one job"
  for job in jobs:
    doAssert job.svtype in ppA.paths, "svtype not in A: " & $job.svtype
    doAssert job.svtype in ppB.paths, "svtype not in B: " & $job.svtype
    doAssert job.chrom in ppA.chromsBySvtype[job.svtype],
      "chrom not in A's set: " & job.chrom & "/" & $job.svtype
    doAssert job.chrom in ppB.chromsBySvtype[job.svtype],
      "chrom not in B's set: " & job.chrom & "/" & $job.svtype

# P06 — job count equals number of (chrom,svtype) pairs in both A and B
timed("P06", "buildWorkQueue: job count equals |intersection of (chrom,svtype)|"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let ppA = preprocessVcf(FixtureA, tmpDir, "A")
  let ppB = preprocessVcf(FixtureB, tmpDir, "B")
  let jobs = buildWorkQueue(ppA, ppB, tmpDir)
  var inBoth = 0
  for svt, chromsA in ppA.chromsBySvtype:
    if svt notin ppB.chromsBySvtype: continue
    for chrom in chromsA:
      if chrom in ppB.chromsBySvtype[svt]:
        inc inBoth
  doAssert jobs.len == inBoth,
    "expected " & $inBoth & " jobs, got " & $jobs.len

# P07 — preprocessVcf works on .bcf inputs
timed("P07", "preprocessVcf: accepts .bcf input"):
  const FixtureA_bcf = "tests/fixtures/fixtureA.bcf"
  doAssert fileExists(FixtureA_bcf), "fixture missing: " & FixtureA_bcf
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA_bcf, tmpDir, "A_bcf")
  doAssert svDEL in pp.paths, "missing DEL from .bcf input"
  doAssert "chr1" in pp.chromsBySvtype[svDEL], "DEL/chr1 should be present"

# Helper: collect every record in a per-SVTYPE BCF as (id, pos, end, svtype, infoNames)
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

# P08 — INFO is slimmed to keep-set only after preprocessing
timed("P08", "preprocessVcf: temp BCF records have only keep-set INFO fields"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  let allowed = ["SVTYPE", "SVLEN", "END", "CHR2", "END2", "POS2",
                 "MATCHA_BOFF"].toHashSet
  for path in pp.paths.values:
    for rec in readRecords(path):
      for name in rec.infoNames:
        doAssert name in allowed,
          "unexpected INFO field " & name & " in " & rec.id

# P09 — record with ID="." gets synthesized (CHROM_POS_SVTYPE_LINENUMBER)
timed("P09", "preprocessVcf: missing ID is synthesized"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  var found = false
  for rec in readRecords(pp.paths[svDEL]):
    # TC14 record sits at chr1:37000 in the input; synthetic ID format is
    # CHROM_POS_SVTYPE_LINENUMBER. The lineno is the 1-based input order.
    if rec.id.startsWith("chr1_37000_DEL_"):
      found = true
      break
  doAssert found, "expected synthesized ID prefix chr1_37000_DEL_<lineno> in DEL BCF"

# P10 — symbolic ALT only (no INFO/SVTYPE) is routed to the right SVTYPE BCF
timed("P10", "preprocessVcf: ALT-symbolic SVTYPE (TC12) routed to DEL BCF"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  var found = false
  for rec in readRecords(pp.paths[svDEL]):
    if rec.id == "DEL_A_12_alt_only":
      doAssert rec.svtype == "DEL", "expected SVTYPE=DEL, got " & rec.svtype
      found = true
      break
  doAssert found, "DEL_A_12_alt_only missing from DEL BCF"

# P11 — SVTYPE conflict (INFO=DUP, ALT=<DEL>) → ALT wins; record lands in DEL BCF
timed("P11", "preprocessVcf: ALT wins on SVTYPE conflict (TC13)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  var inDel = false
  for rec in readRecords(pp.paths[svDEL]):
    if rec.id == "DEL_A_13_conflict":
      doAssert rec.svtype == "DEL", "SVTYPE should be DEL after ALT wins"
      inDel = true
      break
  doAssert inDel, "DEL_A_13_conflict should be in DEL BCF (ALT=<DEL>)"
  if svDUP in pp.paths:
    for rec in readRecords(pp.paths[svDUP]):
      doAssert rec.id != "DEL_A_13_conflict", "should not be in DUP BCF"

# P12 — record with neither END nor SVLEN is skipped
timed("P12", "preprocessVcf: missing END and SVLEN → skipped (TC15)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  for path in pp.paths.values:
    for rec in readRecords(path):
      doAssert rec.id != "DEL_A_15_no_end_no_svlen",
        "TC15 should have been skipped (no END, no SVLEN)"

# P13 — record with END < POS is skipped
timed("P13", "preprocessVcf: END <= POS → skipped (TC16)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  for path in pp.paths.values:
    for rec in readRecords(path):
      doAssert rec.id != "DEL_A_16_bad_end",
        "TC16 should have been skipped (END < POS)"

# P14 — END/SVLEN inconsistency >10% → kept; SVLEN normalized to END-POS
timed("P14", "preprocessVcf: inconsistent END/SVLEN → SVLEN := END-POS (TC17)"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")
  var found = false
  for rec in readRecords(pp.paths[svDEL]):
    if rec.id == "DEL_A_17_inconsistent":
      # Input: POS=43000, END=44000, INFO/SVLEN=-1500. END-POS=1000.
      # |1500-1000|/1500 = 33% > 10% → END wins, SVLEN normalized to 1000.
      doAssert rec.endPos == 44000, "END should be 44000, got " & $rec.endPos
      doAssert rec.svlen == 1000,
        "SVLEN should be normalized to 1000 (END-POS), got " & $rec.svlen
      found = true
      break
  doAssert found, "DEL_A_17_inconsistent missing from DEL BCF"

# P15 — MATCHA_BOFF round-trip: pull the offset from a slim record, seek into
# the original file, decode the record there, confirm CHROM/POS/ID match.
timed("P15", "MATCHA_BOFF: bgzf_seek round-trip into original input"):
  let tmpDir = createTempDir("matcha_test_", "")
  defer: removeDir(tmpDir)
  let pp = preprocessVcf(FixtureA, tmpDir, "A")

  # Find DEL_A_01 in the slim DEL BCF and read its MATCHA_BOFF.
  var slim: VCF
  doAssert open(slim, pp.paths[svDEL]), "cannot open slim DEL BCF"
  var boffData: seq[int32]
  var slimChrom: string
  var slimPos: int64
  var slimId: string
  var found = false
  for v in slim:
    if $v.ID == "DEL_A_01":
      doAssert v.info.get("MATCHA_BOFF", boffData) == Status.OK,
        "MATCHA_BOFF missing on slim DEL_A_01"
      doAssert boffData.len == 2,
        "MATCHA_BOFF should be Number=2, got len=" & $boffData.len
      slimChrom = $v.CHROM
      slimPos = v.POS
      slimId = $v.ID
      found = true
      break
  slim.close()
  doAssert found, "DEL_A_01 missing from slim DEL BCF"

  let offset = (int64(boffData[0]) shl 32) or
               (int64(boffData[1]) and 0xFFFFFFFF'i64)
  doAssert offset > 0, "decoded MATCHA_BOFF should be positive, got " & $offset

  # Seek into the ORIGINAL file at that offset and read one record.
  var orig: VCF
  doAssert open(orig, FixtureA), "cannot open original fixture"
  doAssert bgzf_seek(rawBgzf(orig), offset, 0.cint) >= 0,
    "bgzf_seek failed for offset " & $offset

  var rec = newVariant()
  rec.vcf = orig
  doAssert bcf_read(rawHts(orig), orig.header.hdr, rec.c) >= 0,
    "bcf_read at offset " & $offset & " failed"
  discard bcf_unpack(rec.c, BCF_UN_STR.cint)

  doAssert $rec.CHROM == slimChrom,
    "CHROM mismatch: slim=" & slimChrom & " orig=" & $rec.CHROM
  doAssert rec.POS == slimPos,
    "POS mismatch: slim=" & $slimPos & " orig=" & $rec.POS
  doAssert $rec.ID == slimId,
    "ID mismatch: slim=" & slimId & " orig=" & $rec.ID
  orig.close()
