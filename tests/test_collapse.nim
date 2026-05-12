## Integration tests for matcha collapse.
## Requires the matcha binary (nimble build) and test fixtures (generate_fixtures.py).
## Run from project root: nim c --hints:off -r tests/test_collapse.nim
echo "--------------- Test Collapse ---------------"

import std/[os, osproc, sequtils, strutils]
import test_utils

const BinPath  = "./matcha"
const FixDelly = "tests/fixtures/collapse_delly.vcf.gz"
const FixManta = "tests/fixtures/collapse_manta.vcf.gz"
const FixMulti = "tests/fixtures/collapse_multiallelic.vcf.gz"

proc tmpBcf(): string =
  "/tmp/test_collapse_" & $os.getCurrentProcessId() & ".bcf"

proc run(args: string): (string, int) =
  let t = getEnv("MATCHA_TEST_TIMEOUT", "60")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>/dev/null")

proc runMerged(args: string): (string, int) =
  let t = getEnv("MATCHA_TEST_TIMEOUT", "60")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>&1")

proc viewBcf(path: string): string =
  let (outp, _) = execCmdEx("bcftools view " & path & " 2>/dev/null")
  outp

type ColRecord = object
  chrom, id, filter: string
  pos:               int
  source, sourceList: string
  nSource, nMerged:  int

proc parseCollapsed(vcfText: string): seq[ColRecord] =
  for line in vcfText.splitLines:
    if line.len == 0 or line[0] == '#': continue
    let cols = line.split('\t')
    if cols.len < 8: continue
    var rec = ColRecord(
      chrom:  cols[0],
      pos:    parseInt(cols[1]),
      id:     cols[2],
      filter: cols[6],
    )
    for field in cols[7].split(';'):
      let kv = field.split('=', 1)
      if kv.len != 2: continue
      case kv[0]
      of "SOURCE":     rec.source     = kv[1]
      of "SOURCELIST": rec.sourceList = kv[1]
      of "N_SOURCE":   rec.nSource    = parseInt(kv[1])
      of "N_MERGED":   rec.nMerged    = parseInt(kv[1])
    result.add(rec)

proc collapseRun(extra = ""): (seq[ColRecord], int) =
  ## Run a standard 2-caller collapse → tmpfile, return parsed records + exit code.
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-jaccard 0.5 Delly:" & FixDelly &
                      " Manta:" & FixManta & " -o " & outBcf & " " & extra)
  if code != 0:
    discard tryRemoveFile(outBcf)
    return (@[], code)
  let records = parseCollapsed(viewBcf(outBcf))
  discard tryRemoveFile(outBcf)
  (records, 0)

# ---------------------------------------------------------------------------

# C01 — binary available
timed("C01", "binary available"):
  doAssert fileExists(BinPath), "binary not found: " & BinPath

# C02 — collapse exits 0 and produces a valid BCF
timed("C02", "2-caller collapse: exits 0, BCF readable"):
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-jaccard 0.5 Delly:" & FixDelly &
                      " Manta:" & FixManta & " -o " & outBcf)
  doAssert code == 0, "collapse exited " & $code
  doAssert fileExists(outBcf), "output BCF not created"
  let view = viewBcf(outBcf)
  doAssert "SVTYPE" in view, "BCF appears empty or malformed"
  discard tryRemoveFile(outBcf)

# C03 — correct record count: 5 output records for 2-caller fixture
timed("C03", "correct output record count (5 for the 2-caller fixture)"):
  let (recs, code) = collapseRun()
  doAssert code == 0, "collapse failed"
  doAssert recs.len == 5, "expected 5 records, got " & $recs.len &
    ": " & recs.mapIt(it.id).join(", ")

# C04 — PASS filter beats LowQual: DEL_M_02 chosen over DEL_D_02
timed("C04", "PASS priority: Manta DEL_M_02 preferred over LowQual DEL_D_02"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  var found = false
  for r in recs:
    if r.pos == 3100 and r.source == "Manta":
      doAssert r.id == "DEL_M_02", "expected DEL_M_02, got " & r.id
      doAssert r.filter == "PASS", "expected PASS, got " & r.filter
      found = true
  doAssert found, "no PASS-selected Manta record at pos 3100"

# C05 — ORDER tiebreak when both PASS: Delly (callerIdx 0) wins over Manta
timed("C05", "ORDER tiebreak: Delly DEL_D_01 preferred over Manta DEL_M_01"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  var found = false
  for r in recs:
    if r.pos == 1000 and r.nMerged == 2:
      doAssert r.source == "Delly", "expected source=Delly, got " & r.source
      doAssert r.id == "DEL_D_01", "expected DEL_D_01, got " & r.id
      found = true
  doAssert found, "no merged record at pos 1000"

# C06 — singletons preserved: DEL_D_03 and DEL_M_03 both appear
timed("C06", "singletons preserved with N_MERGED=1 and N_SOURCE=1"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  var d03, m03: bool
  for r in recs:
    if r.id == "DEL_D_03":
      doAssert r.nMerged == 1, "DEL_D_03 nMerged=" & $r.nMerged
      doAssert r.nSource == 1, "DEL_D_03 nSource=" & $r.nSource
      doAssert r.sourceList == "Delly", "DEL_D_03 sourceList=" & r.sourceList
      d03 = true
    elif r.id == "DEL_M_03":
      doAssert r.nMerged == 1, "DEL_M_03 nMerged=" & $r.nMerged
      doAssert r.nSource == 1, "DEL_M_03 nSource=" & $r.nSource
      doAssert r.sourceList == "Manta", "DEL_M_03 sourceList=" & r.sourceList
      m03 = true
  doAssert d03, "DEL_D_03 singleton not found"
  doAssert m03, "DEL_M_03 singleton not found"

# C07 — merged records carry N_MERGED=2 and N_SOURCE=2
timed("C07", "merged records: N_MERGED=2, N_SOURCE=2, SOURCELIST contains both"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  for r in recs:
    if r.nMerged == 2:
      doAssert r.nSource == 2, r.id & " nSource=" & $r.nSource
      doAssert "Delly" in r.sourceList and "Manta" in r.sourceList,
               r.id & " sourceList=" & r.sourceList

# C08 — DUP cluster: Delly wins (ORDER), N_MERGED=2
timed("C08", "DUP cluster: DUP_D_01 selected by ORDER, N_MERGED=2"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  var found = false
  for r in recs:
    if r.id == "DUP_D_01":
      doAssert r.source == "Delly"
      doAssert r.nMerged == 2
      found = true
  doAssert found, "DUP_D_01 not found in output"

# C09 — chromosomal sort order: records sorted by pos within chrom
timed("C09", "output records sorted by position"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  for i in 1 ..< recs.len:
    doAssert recs[i].pos >= recs[i-1].pos,
      "out of order: " & $recs[i-1].pos & " then " & $recs[i].pos

# C10 — multiallelic input causes error exit
timed("C10", "multiallelic record causes non-zero exit with informative message"):
  let outBcf = tmpBcf()
  let (outp, code) = runMerged("collapse --min-jaccard 0.5 Caller:" & FixMulti &
                                " -o " & outBcf)
  discard tryRemoveFile(outBcf)
  doAssert code != 0, "expected non-zero exit for multiallelic input, got 0"
  doAssert "multiallelic" in outp.toLowerAscii,
    "expected 'multiallelic' in error output, got: " & outp

# C11 — --min-overlap threshold also works (not only --min-jaccard)
timed("C11", "--min-overlap 0.5 produces the same cluster topology"):
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-overlap 0.5 Delly:" & FixDelly &
                      " Manta:" & FixManta & " -o " & outBcf)
  doAssert code == 0, "collapse --min-overlap failed"
  let recs = parseCollapsed(viewBcf(outBcf))
  discard tryRemoveFile(outBcf)
  doAssert recs.len == 5, "expected 5 records with --min-overlap 0.5, got " & $recs.len

# C12 — --info filter: requesting only SVTYPE in output drops END/SVLEN
timed("C12", "--info SVTYPE filter keeps SVTYPE, drops END and SVLEN"):
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-jaccard 0.5 --info SVTYPE Delly:" & FixDelly &
                      " Manta:" & FixManta & " -o " & outBcf)
  doAssert code == 0, "collapse --info SVTYPE failed"
  let view = viewBcf(outBcf)
  discard tryRemoveFile(outBcf)
  var sawSvtype, sawEnd, sawSvlen: bool
  for line in view.splitLines:
    if line.len == 0 or line[0] == '#': continue
    if "SVTYPE=" in line: sawSvtype = true
    if "END=" in line:    sawEnd    = true
    if "SVLEN=" in line:  sawSvlen  = true
  doAssert sawSvtype, "--info SVTYPE: SVTYPE not in output"
  doAssert not sawEnd,   "--info SVTYPE: END should be absent, but found"
  doAssert not sawSvlen, "--info SVTYPE: SVLEN should be absent, but found"
