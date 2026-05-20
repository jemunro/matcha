## Integration tests for matcha merge.
## Requires the matcha binary (nimble build) and test fixtures
## (generate_fixtures.py). Run from project root:
##   nim c --hints:off -r tests/test_merge.nim
echo "--------------- Test Merge ---------------"

import std/[os, osproc, sequtils, strutils, tables]
import test_utils

const BinPath = "./matcha"
const FixS1     = "tests/fixtures/merge_S1.vcf.gz"
const FixS2     = "tests/fixtures/merge_S2.vcf.gz"
const FixS3     = "tests/fixtures/merge_S3.vcf.gz"
const FixS1Dup  = "tests/fixtures/merge_S1_dup.vcf.gz"
const Fix2Samp  = "tests/fixtures/merge_2sample.vcf.gz"

proc tmpOut(ext = ".vcf"): string =
  "/tmp/test_merge_" & $os.getCurrentProcessId() & ext

proc run(args: string): (string, int) =
  let t = getEnv("MATCHA_TEST_TIMEOUT", "60")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>/dev/null")

proc runMerged(args: string): (string, int) =
  let t = getEnv("MATCHA_TEST_TIMEOUT", "60")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>&1")

type MergeRec = object
  chrom, id, alt:    string
  pos:               int
  ac, an:            int
  af:                float
  gts:               seq[string]   ## per-sample GT in CLI sample-column order

proc parseInfoKV(info: string): Table[string, string] =
  for f in info.split(';'):
    let kv = f.split('=', 1)
    if kv.len == 2: result[kv[0]] = kv[1]

proc parseMergeVcf(text: string): tuple[samples: seq[string]; rows: seq[MergeRec]] =
  for line in text.splitLines:
    if line.len == 0: continue
    if line.startsWith("#CHROM"):
      let cols = line.split('\t')
      if cols.len >= 10:
        result.samples = cols[9 .. ^1]
      continue
    if line[0] == '#': continue
    let cols = line.split('\t')
    if cols.len < 9: continue
    var r = MergeRec(
      chrom: cols[0], pos: parseInt(cols[1]),
      id: cols[2], alt: cols[4],
    )
    let info = parseInfoKV(cols[7])
    if "AC" in info: r.ac = parseInt(info["AC"])
    if "AN" in info: r.an = parseInt(info["AN"])
    if "AF" in info:
      try: r.af = parseFloat(info["AF"]) except ValueError: r.af = 0.0
    # FORMAT in col 8; per-sample GT in cols 9..
    let fmtTokens = cols[8].split(':')
    var gtIdx = -1
    for i, t in fmtTokens.pairs:
      if t == "GT": gtIdx = i; break
    for s in cols[9 .. ^1]:
      let stoks = s.split(':')
      if gtIdx >= 0 and gtIdx < stoks.len:
        r.gts.add(stoks[gtIdx])
      else:
        r.gts.add(".")
    result.rows.add(r)

proc mergeRun(extra = ""; outExt = ".vcf"): (seq[string], seq[MergeRec], int) =
  let outPath = tmpOut(outExt)
  let (_, code) = run("merge --min-jaccard 0.75 " &
                      "S1:" & FixS1 & " S2:" & FixS2 & " S3:" & FixS3 &
                      " -o " & outPath & " " & extra)
  if code != 0:
    discard tryRemoveFile(outPath)
    return (@[], @[], code)
  let text = readFile(outPath)
  discard tryRemoveFile(outPath)
  let parsed = parseMergeVcf(text)
  (parsed.samples, parsed.rows, 0)

# ---------------------------------------------------------------------------

# M01 — binary available
timed("M01", "binary available"):
  doAssert fileExists(BinPath), "binary not found: " & BinPath

# M02 — merge exits 0 and produces output with 3 samples
timed("M02", "3-sample merge: exits 0 with 3 sample columns"):
  let (samples, rows, code) = mergeRun()
  doAssert code == 0, "merge exited " & $code
  doAssert samples == @["S1", "S2", "S3"],
           "expected sample columns S1,S2,S3 in order, got " & samples.join(",")
  doAssert rows.len == 6, "expected 6 cluster rows, got " & $rows.len

# M03 — DEL_1000 cluster: all 3 samples, AC=3 AN=6 AF=0.5
timed("M03", "DEL_1000 cluster: AC=3 AN=6 AF=0.5 with S1=0/1 S2=1/1 S3=0/0"):
  let (_, rows, code) = mergeRun()
  doAssert code == 0
  var found = false
  for r in rows:
    if r.pos == 1000:
      doAssert r.ac == 3, "AC=" & $r.ac
      doAssert r.an == 6, "AN=" & $r.an
      doAssert abs(r.af - 0.5) < 1e-3, "AF=" & $r.af
      doAssert r.gts == @["0/1", "1/1", "0/0"], "GTs=" & r.gts.join(",")
      found = true
  doAssert found, "no row at pos 1000"

# M04 — DEL_3000 cluster: only S1, others ./. ; AC=2 AN=2 AF=1.0
timed("M04", "DEL_3000 singleton: only S1 contributes; AN=2 AC=2 AF=1.0"):
  let (_, rows, code) = mergeRun()
  doAssert code == 0
  var found = false
  for r in rows:
    if r.pos == 3000:
      doAssert r.ac == 2, "AC=" & $r.ac
      doAssert r.an == 2, "AN=" & $r.an
      doAssert abs(r.af - 1.0) < 1e-3, "AF=" & $r.af
      doAssert r.gts == @["1/1", "./.", "./."], "GTs=" & r.gts.join(",")
      found = true
  doAssert found, "no row at pos 3000"

# M05 — DUP_9000 cluster: S1+S2, S3 missing; AC=1 AN=4 AF=0.25
timed("M05", "DUP_9000: S1=0/1 S2=0/0 S3=./. ; AC=1 AN=4"):
  let (_, rows, code) = mergeRun()
  doAssert code == 0
  for r in rows:
    if r.pos == 9000:
      doAssert r.ac == 1
      doAssert r.an == 4
      doAssert r.gts == @["0/1", "0/0", "./."], r.gts.join(",")

# M06 — INV_28000: only S3; AC=1 AN=2 AF=0.5
timed("M06", "INV_28000: only S3 contributes; AN=2 AC=1"):
  let (_, rows, code) = mergeRun()
  doAssert code == 0
  for r in rows:
    if r.pos == 28000:
      doAssert r.ac == 1
      doAssert r.an == 2
      doAssert r.gts == @["./.", "./.", "0/1"], r.gts.join(",")

# M07 — BND cluster: all 3 samples; AC=4 AN=6 AF=2/3; ALT preserved
timed("M07", "BND cluster: ALT bracket-form preserved; AC=4 AN=6"):
  let (_, rows, code) = mergeRun()
  doAssert code == 0
  var found = false
  for r in rows:
    if r.pos == 31000:
      doAssert r.ac == 4, "AC=" & $r.ac
      doAssert r.an == 6
      doAssert "[chr1:32000[" in r.alt or "]chr1:32000]" in r.alt or
               "[chr1:32000" in r.alt,
               "BND ALT not preserved, got: " & r.alt
      doAssert r.gts == @["0/1", "0/1", "1/1"], r.gts.join(",")
      found = true
  doAssert found, "no BND row at pos 31000"

# M08 — multi-sample input rejected
timed("M08", "input with >1 sample columns is rejected"):
  let outPath = tmpOut()
  let (outp, code) = runMerged("merge " & Fix2Samp & " " & FixS2 & " -o " & outPath)
  discard tryRemoveFile(outPath)
  doAssert code != 0, "expected non-zero exit for multi-sample input"
  doAssert "sample column" in outp.toLowerAscii or "exactly 1" in outp,
           "expected sample-count error, got: " & outp

# M09 — duplicate sample IDs rejected
timed("M09", "duplicate sample IDs across inputs rejected"):
  let outPath = tmpOut()
  let (outp, code) = runMerged("merge " & FixS1 & " " & FixS1Dup & " -o " & outPath)
  discard tryRemoveFile(outPath)
  doAssert code != 0, "expected non-zero exit for duplicate sample IDs"
  doAssert "duplicate sample" in outp.toLowerAscii,
           "expected duplicate-sample-id error, got: " & outp

# M10 — --format omitting GT still produces cohort AC/AN/AF (auto-add)
timed("M10", "--format DP (no GT) silently adds GT; cohort fields still correct"):
  let outPath = tmpOut()
  let (_, code) = run("merge --min-jaccard 0.75 --format DP " &
                      "S1:" & FixS1 & " S2:" & FixS2 & " S3:" & FixS3 &
                      " -o " & outPath)
  doAssert code == 0, "merge with --format DP failed"
  let text = readFile(outPath)
  discard tryRemoveFile(outPath)
  let parsed = parseMergeVcf(text)
  doAssert parsed.rows.len == 6
  for r in parsed.rows:
    if r.pos == 1000:
      doAssert r.ac == 3
      doAssert r.an == 6

# M11 — sort order: rows in (chrom, pos) order
timed("M11", "output rows sorted by (chrom, pos)"):
  let (_, rows, code) = mergeRun()
  doAssert code == 0
  for i in 1 ..< rows.len:
    doAssert rows[i].pos >= rows[i-1].pos,
             "out of order: " & $rows[i-1].pos & " then " & $rows[i].pos

# M12 — output INFO does not leak matcha-internal fields
timed("M12", "no SRC_INDEX, CALLER_IDX, or FORMAT/SID in output"):
  let outPath = tmpOut()
  let (_, code) = run("merge --min-jaccard 0.75 " &
                      "S1:" & FixS1 & " S2:" & FixS2 & " S3:" & FixS3 &
                      " -o " & outPath)
  doAssert code == 0
  let text = readFile(outPath)
  discard tryRemoveFile(outPath)
  doAssert "SRC_INDEX" notin text, "SRC_INDEX leaked into output"
  doAssert "CALLER_IDX" notin text, "CALLER_IDX leaked into output"
  doAssert "SID" notin text, "FORMAT/SID leaked into output"

# M13 — < 2 inputs rejected
timed("M13", "single input file rejected"):
  let outPath = tmpOut()
  let (outp, code) = runMerged("merge " & FixS1 & " -o " & outPath)
  discard tryRemoveFile(outPath)
  doAssert code != 0
  doAssert "at least 2" in outp.toLowerAscii or "got 1" in outp,
           "expected min-inputs error, got: " & outp

# M14 — N_MERGED is not present (cohort context)
timed("M14", "cohort output does not emit N_MERGED"):
  let outPath = tmpOut()
  let (_, code) = run("merge --min-jaccard 0.75 " &
                      "S1:" & FixS1 & " S2:" & FixS2 & " S3:" & FixS3 &
                      " -o " & outPath)
  doAssert code == 0
  let text = readFile(outPath)
  discard tryRemoveFile(outPath)
  doAssert "N_MERGED" notin text, "N_MERGED appeared in output"

# M15 — INS cluster: S1 and S2 contribute (S1=0/1, S2=1/1, S3 missing).
# Sequence-resolved ALT is preserved on the merged record; AC=3, AN=4, S3=./..
timed("M15", "INS cluster: sequence ALT preserved, AC=3 AN=4 with S1/S2 only"):
  let (samples, rows, code) = mergeRun()
  doAssert code == 0
  var found = false
  for r in rows:
    if r.chrom == "chr1" and r.pos == 40000:
      found = true
      doAssert r.alt.startsWith("N") and r.alt.len >= 50,
        "INS ALT should be sequence-resolved, got: '" & r.alt & "'"
      doAssert "<INS>" notin r.alt,
        "INS ALT should not have collapsed to <INS>, got: '" & r.alt & "'"
      doAssert r.ac == 3, "INS AC expected 3, got " & $r.ac
      doAssert r.an == 4, "INS AN expected 4, got " & $r.an
  doAssert found, "INS cluster at chr1:40000 missing from merge output"
