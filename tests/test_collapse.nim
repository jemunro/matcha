## Integration tests for matcha collapse.
## Requires the matcha binary (nimble build) and test fixtures (generate_fixtures.py).
## Run from project root: nim c --hints:off -r tests/test_collapse.nim
echo "--------------- Test Collapse ---------------"

import std/[os, osproc, sequtils, strutils, tables]
import test_utils

const BinPath    = "./matcha"
const FixCaller1   = "tests/fixtures/collapse_caller1.vcf.gz"
const FixCaller2   = "tests/fixtures/collapse_caller2.vcf.gz"
const FixMulti   = "tests/fixtures/collapse_multiallelic.vcf.gz"
const FixCaller1_1S = "tests/fixtures/collapse_caller1_1sample.vcf.gz"
const FixCaller2_1S = "tests/fixtures/collapse_caller2_1sample.vcf.gz"
const FixCaller1_2S = "tests/fixtures/collapse_caller1_2sample.vcf.gz"

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
  let (_, code) = run("collapse --min-jaccard 0.5 Caller1:" & FixCaller1 &
                      " Caller2:" & FixCaller2 & " -o " & outBcf & " " & extra)
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
  let (_, code) = run("collapse --min-jaccard 0.5 Caller1:" & FixCaller1 &
                      " Caller2:" & FixCaller2 & " -o " & outBcf)
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
timed("C04", "PASS priority: Caller2 DEL_M_02 preferred over LowQual DEL_D_02"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  var found = false
  for r in recs:
    if r.pos == 3100 and r.source == "Caller2":
      doAssert r.id == "DEL_M_02", "expected DEL_M_02, got " & r.id
      doAssert r.filter == "PASS", "expected PASS, got " & r.filter
      found = true
  doAssert found, "no PASS-selected Caller2 record at pos 3100"

# C05 — ORDER tiebreak when both PASS: Caller1 (callerIdx 0) wins over Caller2
timed("C05", "ORDER tiebreak: Caller1 DEL_D_01 preferred over Caller2 DEL_M_01"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  var found = false
  for r in recs:
    if r.pos == 1000 and r.nMerged == 2:
      doAssert r.source == "Caller1", "expected source=Caller1, got " & r.source
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
      doAssert r.sourceList == "Caller1", "DEL_D_03 sourceList=" & r.sourceList
      d03 = true
    elif r.id == "DEL_M_03":
      doAssert r.nMerged == 1, "DEL_M_03 nMerged=" & $r.nMerged
      doAssert r.nSource == 1, "DEL_M_03 nSource=" & $r.nSource
      doAssert r.sourceList == "Caller2", "DEL_M_03 sourceList=" & r.sourceList
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
      doAssert "Caller1" in r.sourceList and "Caller2" in r.sourceList,
               r.id & " sourceList=" & r.sourceList

# C08 — DUP cluster: Caller1 wins (ORDER), N_MERGED=2
timed("C08", "DUP cluster: DUP_D_01 selected by ORDER, N_MERGED=2"):
  let (recs, code) = collapseRun()
  doAssert code == 0
  var found = false
  for r in recs:
    if r.id == "DUP_D_01":
      doAssert r.source == "Caller1"
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
  let (_, code) = run("collapse --min-overlap 0.5 Caller1:" & FixCaller1 &
                      " Caller2:" & FixCaller2 & " -o " & outBcf)
  doAssert code == 0, "collapse --min-overlap failed"
  let recs = parseCollapsed(viewBcf(outBcf))
  discard tryRemoveFile(outBcf)
  doAssert recs.len == 5, "expected 5 records with --min-overlap 0.5, got " & $recs.len

# C12 — --info filter: requesting only SVTYPE in output drops END/SVLEN
timed("C12", "--info SVTYPE filter keeps SVTYPE, drops END and SVLEN"):
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-jaccard 0.5 --info SVTYPE Caller1:" & FixCaller1 &
                      " Caller2:" & FixCaller2 & " -o " & outBcf)
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

# CT03 — --format "" → output has zero sample columns
timed("CT03", "--format '' produces zero-sample output"):
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-jaccard 0.5 --format \"\" Caller1:" & FixCaller1 &
                      " Caller2:" & FixCaller2 & " -o " & outBcf)
  doAssert code == 0, "collapse --format '' failed"
  let view = viewBcf(outBcf)
  discard tryRemoveFile(outBcf)
  # #CHROM header line should have exactly 8 columns (no FORMAT, no samples).
  var sawChromLine = false
  for line in view.splitLines:
    if line.startsWith("#CHROM"):
      sawChromLine = true
      let n = line.split('\t').len
      doAssert n == 8, "expected 8 columns in #CHROM line (no FORMAT/samples), got " & $n
  doAssert sawChromLine, "no #CHROM line in output"

# CT05 — 1-sample collapse: output has one sample column (caller 0's sample
# name) and GT carries through per record from its source caller.
timed("CT05", "1-sample collapse: SAMPLE1 column carries GT round-trip"):
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-jaccard 0.5 Caller1:" & FixCaller1_1S &
                      " Caller2:" & FixCaller2_1S & " -o " & outBcf)
  doAssert code == 0, "1-sample collapse failed"
  let view = viewBcf(outBcf)
  discard tryRemoveFile(outBcf)
  # #CHROM line should have 10 cols (8 std + FORMAT + 1 sample); sample is SAMPLE1.
  var chromCols: seq[string]
  for line in view.splitLines:
    if line.startsWith("#CHROM"):
      chromCols = line.split('\t')
      break
  doAssert chromCols.len == 10, "expected 10 cols in #CHROM, got " & $chromCols.len
  doAssert chromCols[8] == "FORMAT", "col 8 = FORMAT, got " & chromCols[8]
  doAssert chromCols[9] == "SAMPLE1", "col 9 = SAMPLE1, got " & chromCols[9]
  # Verify GTs round-trip. Per fixture: DEL_D_01 (Caller1)=0/1, DEL_M_02 (Caller2)=0/1,
  # DEL_D_03 (Caller1)=1/1, DEL_M_03 (Caller2)=0/0, DUP_D_01 (Caller1)=0/1.
  var gts: Table[string, string]
  for line in view.splitLines:
    if line.len == 0 or line[0] == '#': continue
    let cols = line.split('\t')
    if cols.len < 10: continue
    gts[cols[2]] = cols[9]
  doAssert gts.getOrDefault("DEL_D_01") == "0/1",
           "DEL_D_01 GT=" & gts.getOrDefault("DEL_D_01")
  doAssert gts.getOrDefault("DEL_M_02") == "0/1",
           "DEL_M_02 GT=" & gts.getOrDefault("DEL_M_02")
  doAssert gts.getOrDefault("DEL_D_03") == "1/1",
           "DEL_D_03 GT=" & gts.getOrDefault("DEL_D_03")
  doAssert gts.getOrDefault("DEL_M_03") == "0/0",
           "DEL_M_03 GT=" & gts.getOrDefault("DEL_M_03")
  doAssert gts.getOrDefault("DUP_D_01") == "0/1",
           "DUP_D_01 GT=" & gts.getOrDefault("DUP_D_01")

# CT06 — 2-sample input rejected with a clear error.
timed("CT06", "multi-sample input causes non-zero exit with informative message"):
  let outBcf = tmpBcf()
  let (outp, code) = runMerged("collapse --min-jaccard 0.5 D:" & FixCaller1_2S &
                               " M:" & FixCaller2_1S & " -o " & outBcf)
  discard tryRemoveFile(outBcf)
  doAssert code != 0, "expected non-zero exit for 2-sample input, got 0"
  doAssert "sample" in outp.toLowerAscii,
           "expected 'sample' in error output, got: " & outp
  doAssert "at most 1" in outp or "supports" in outp,
           "expected informative reject message, got: " & outp

# CT07 — inconsistent sample counts across inputs rejected.
timed("CT07", "inconsistent sample counts cause non-zero exit"):
  let outBcf = tmpBcf()
  let (outp, code) = runMerged("collapse --min-jaccard 0.5 D:" & FixCaller1_1S &
                               " M:" & FixCaller2 & " -o " & outBcf)
  discard tryRemoveFile(outBcf)
  doAssert code != 0, "expected non-zero exit for inconsistent counts, got 0"
  doAssert "inconsistent" in outp.toLowerAscii or "same sample count" in outp,
           "expected inconsistent-counts message, got: " & outp

# CT04 — 3-caller run (aliasing one fixture twice) — streaming loop ticks
# across 3 readers; record counts/sources are correct.
timed("CT04", "3-caller run: streaming across 3 readers produces sensible counts"):
  let outBcf = tmpBcf()
  let (_, code) = run("collapse --min-jaccard 0.5 " &
                      " Caller1:" & FixCaller1 &
                      " Caller2:" & FixCaller2 &
                      " Caller1b:" & FixCaller1 &
                      " -o " & outBcf)
  doAssert code == 0, "3-caller collapse failed"
  let recs = parseCollapsed(viewBcf(outBcf))
  discard tryRemoveFile(outBcf)
  # The Caller1 fixture self-matches Caller1b (identical content) → most records
  # cluster across Caller1 and Caller1b with N_SOURCE >= 2.
  doAssert recs.len > 0, "no records in 3-caller output"
  var sawCaller1b: bool
  for r in recs:
    if "Caller1b" in r.sourceList: sawCaller1b = true
  doAssert sawCaller1b, "Caller1b (third caller) absent from any cluster SOURCELIST"
