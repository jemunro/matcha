## Tests for matcha anno — expression parsing, aggregation kernels, and
## end-to-end annotation against the DB fixture.
## Run from project root: nim c --hints:off -r tests/test_anno.nim
echo "--------------- Test Anno ---------------"

import std/[os, osproc, strutils, tables]
import test_utils
import matcha/anno
import matcha/log

setQuiet(true)

const BinPath   = "./matcha"
const FixtureA  = "tests/fixtures/fixtureA.vcf.gz"
const FixtureDB = "tests/fixtures/fixtureDB.vcf.gz"

proc run(args: string): (string, int) =
  ## Memoized — identical commands within a single test process share output.
  runMatcha(BinPath, args, getEnv("MATCHA_TEST_TIMEOUT", "30"))

proc runMerged(args: string): (string, int) =
  let t = getEnv("MATCHA_TEST_TIMEOUT", "30")
  execCmdEx("timeout " & t & " " & BinPath & " " & args & " 2>&1")

# Parse a VCF text into (id → info-dict). Only INFO key=value pairs make it
# in; flags (FOO without =) are stored as "FOO" → "".
proc parseInfo(vcfText: string): Table[string, Table[string, string]] =
  for line in vcfText.strip.splitLines:
    if line.len == 0 or line[0] == '#': continue
    let cols = line.split('\t')
    if cols.len < 8: continue
    let id = cols[2]
    var fields: Table[string, string]
    for kv in cols[7].split(';'):
      let eq = kv.find('=')
      if eq < 0: fields[kv] = ""
      else: fields[kv[0 ..< eq]] = kv[eq + 1 .. ^1]
    result[id] = fields

# ---------------------------------------------------------------------------
# Expression parsing (unit, no I/O)
# ---------------------------------------------------------------------------

timed("A01", "parseAnnoExpr: OUTFIELD=FUNC(SRCFIELD) round-trip"):
  let e = parseAnnoExpr("AF_DB=max(AF)")
  doAssert e.outField == "AF_DB"
  doAssert e.fn == afMax
  doAssert e.srcField == "AF"
  doAssert e.matchaVar == mvNone

timed("A02", "parseAnnoExpr: unknown FUNC raises ValueError"):
  var raised = false
  try: discard parseAnnoExpr("X=median(Y)")
  except ValueError: raised = true
  doAssert raised, "median should not be a valid FUNC"

timed("A03", "parseAnnoExpr: malformed input raises"):
  for bad in [ "no_equals(x)", "BAD=", "BAD=funcwithoutparens",
               "BAD=max()", "=max(X)", "BAD=max(BAD NAME)" ]:
    var raised = false
    try: discard parseAnnoExpr(bad)
    except ValueError: raised = true
    doAssert raised, "expected parse error for: " & bad

timed("A04", "parseAnnoExpr: MATCHA_* SRCFIELDs set matchaVar correctly"):
  doAssert parseAnnoExpr("N=first(MATCHA_COUNT)").matchaVar == mvCount
  doAssert parseAnnoExpr("S=max(MATCHA_SIMILARITY)").matchaVar == mvSimilarity
  # The old MATCHA_OVERLAP / MATCHA_JACCARD names are no longer recognized;
  # they fall through to mvNone and then fail validation as unknown DB fields.
  doAssert parseAnnoExpr("J=max(MATCHA_JACCARD)").matchaVar == mvNone

# ---------------------------------------------------------------------------
# Aggregation kernel (unit, hand-crafted match sets)
# ---------------------------------------------------------------------------

proc mkMatch(posB: int64; sim: float64;
             payload: Table[string, seq[string]] = initTable[string, seq[string]]()
            ): AnnoMatch =
  AnnoMatch(aIndex: 0, posB: posB,
            similarity: sim, payload: payload)

timed("A05", "applyAggFunc: max/min/mean over numeric DB field"):
  var e = parseAnnoExpr("X=max(AF)"); e.dbType = "Float"; e.dbNumber = "1"
  let m1 = mkMatch(100, 1.0, {"AF": @["0.1"]}.toTable)
  let m2 = mkMatch(200, 1.0, {"AF": @["0.5"]}.toTable)
  let m3 = mkMatch(300, 1.0, {"AF": @["0.3"]}.toTable)
  let matches = @[m1, m2, m3]
  doAssert applyAggFunc(e, matches)[0].startsWith("0.5")

  e.fn = afMin
  doAssert applyAggFunc(e, matches)[0].startsWith("0.1")

  e.fn = afMean
  let mean = parseFloat(applyAggFunc(e, matches)[0])
  doAssert abs(mean - 0.3) < 1e-5, "mean should be 0.3, got " & $mean

timed("A06", "applyAggFunc: first/last respect posB ordering"):
  var e = parseAnnoExpr("X=first(CALLER)"); e.dbType = "String"; e.dbNumber = "1"
  let matches = @[
    mkMatch(300, 1.0, {"CALLER": @["c"]}.toTable),
    mkMatch(100, 1.0, {"CALLER": @["a"]}.toTable),
    mkMatch(200, 1.0, {"CALLER": @["b"]}.toTable),
  ]
  doAssert applyAggFunc(e, matches) == @["a"]
  e.fn = afLast
  doAssert applyAggFunc(e, matches) == @["c"]

timed("A07", "applyAggFunc: best ranks by similarity"):
  var e = parseAnnoExpr("X=best(CALLER)"); e.dbType = "String"; e.dbNumber = "1"
  let matches = @[
    mkMatch(100, 0.5, {"CALLER": @["lo"]}.toTable),
    mkMatch(200, 0.95, {"CALLER": @["top"]}.toTable),
    mkMatch(300, 0.6, {"CALLER": @["mid"]}.toTable),
  ]
  doAssert applyAggFunc(e, matches) == @["top"]

timed("A08", "applyAggFunc: best tie-break by earliest posB"):
  var e = parseAnnoExpr("X=best(CALLER)"); e.dbType = "String"; e.dbNumber = "1"
  let matches = @[
    mkMatch(300, 0.9, {"CALLER": @["late"]}.toTable),
    mkMatch(100, 0.9, {"CALLER": @["early"]}.toTable),
    mkMatch(200, 0.9, {"CALLER": @["mid"]}.toTable),
  ]
  doAssert applyAggFunc(e, matches) == @["early"]

timed("A09", "applyAggFunc: all flattens list-valued source across matches"):
  var e = parseAnnoExpr("X=all(CALLERS)"); e.dbType = "String"; e.dbNumber = "."
  let matches = @[
    mkMatch(100, 1.0, {"CALLERS": @["a", "b"]}.toTable),
    mkMatch(200, 1.0, {"CALLERS": @["c"]}.toTable),
  ]
  doAssert applyAggFunc(e, matches) == @["a", "b", "c"]

timed("A10", "applyAggFunc: unique deduplicates flattened list"):
  var e = parseAnnoExpr("X=unique(CALLERS)"); e.dbType = "String"; e.dbNumber = "."
  let matches = @[
    mkMatch(100, 1.0, {"CALLERS": @["a", "b"]}.toTable),
    mkMatch(200, 1.0, {"CALLERS": @["b", "c"]}.toTable),
  ]
  doAssert applyAggFunc(e, matches) == @["a", "b", "c"]

timed("A11", "applyAggFunc: empty match set returns absent except MATCHA_COUNT=0"):
  var e1 = parseAnnoExpr("X=max(AF)"); e1.dbType = "Float"; e1.dbNumber = "1"
  doAssert applyAggFunc(e1, @[]).len == 0

  let e2 = parseAnnoExpr("N=first(MATCHA_COUNT)")
  doAssert applyAggFunc(e2, @[]) == @["0"]

  let e3 = parseAnnoExpr("S=max(MATCHA_SIMILARITY)")
  doAssert applyAggFunc(e3, @[]).len == 0

timed("A12", "applyAggFunc: MATCHA_SIMILARITY vector vs scalar MATCHA_COUNT"):
  let matches = @[
    mkMatch(100, 0.9, initTable[string, seq[string]]()),
    mkMatch(200, 0.4, initTable[string, seq[string]]()),
  ]
  let eMax = parseAnnoExpr("S=max(MATCHA_SIMILARITY)")
  doAssert parseFloat(applyAggFunc(eMax, matches)[0]) - 0.9 < 1e-5

  let eCount = parseAnnoExpr("N=first(MATCHA_COUNT)")
  doAssert applyAggFunc(eCount, matches) == @["2"]

# ---------------------------------------------------------------------------
# End-to-end integration tests against the matcha binary
# ---------------------------------------------------------------------------

timed("A13", "binary available with anno subcommand"):
  if not fileExists(BinPath):
    let (outp, code) = execCmdEx("nimble build 2>&1")
    if code != 0:
      echo "nimble build failed:\n", outp; quit(1)
  let (outp, code) = runMerged("anno --help")
  doAssert code == 0, "anno --help exit " & $code
  doAssert "OUTFIELD=FUNC(SRCFIELD)" in outp

timed("A14", "no -a expression: exits non-zero"):
  let (outp, code) = runMerged("anno --min-overlap 0.5 " & FixtureA & " " & FixtureDB)
  doAssert code != 0
  doAssert "-a" in outp or "expression" in outp,
    "error should mention missing -a: " & outp

timed("A15", "no threshold: uses default --min-jaccard 0.75"):
  let (outp, code) = run("anno -a X=max\\(AF\\) " & FixtureA & " " & FixtureDB)
  doAssert code == 0, "should exit 0 with default threshold, got: " & outp

timed("A16", "unknown DB SRCFIELD: exits non-zero"):
  let (outp, code) = runMerged("anno -a X=max\\(NOPE\\) --min-overlap 0.5 " &
                                FixtureA & " " & FixtureDB)
  doAssert code != 0
  doAssert "NOPE" in outp

timed("A17", "end-to-end: max/min/mean/unique annotations"):
  let (outp, code) = run("anno --min-overlap 0.5 " &
    "-a AF_MAX=max\\(AF\\) -a AF_MIN=min\\(AF\\) -a AF_MEAN=mean\\(AF\\) " &
    "-a CALLERS_U=unique\\(CALLERS\\) -a N=first\\(MATCHA_COUNT\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code == 0, "exit " & $code
  let info = parseInfo(outp)
  doAssert "DEL_A_06" in info
  let a06 = info["DEL_A_06"]
  # DEL_A_06 matches DEL_DB_06a (AF=0.2) and DEL_DB_06b (AF=0.3)
  doAssert a06["AF_MAX"].startsWith("0.3"), "AF_MAX=" & a06["AF_MAX"]
  doAssert a06["AF_MIN"].startsWith("0.2"), "AF_MIN=" & a06["AF_MIN"]
  doAssert parseFloat(a06["AF_MEAN"]) - 0.25 < 1e-5, "AF_MEAN=" & a06["AF_MEAN"]
  # CALLERS: DB 06b is at posB=12950 (sorts first), DB 06a at 13050.
  # Flattened: [delly,gatk]+[manta] → unique: delly,gatk,manta.
  doAssert a06["CALLERS_U"] == "delly,gatk,manta", "got " & a06["CALLERS_U"]
  doAssert a06["N"] == "2", "N=" & a06["N"]

timed("A18", "end-to-end: MATCHA_COUNT=0 on unmatched, absent fields"):
  let (outp, code) = run("anno --min-overlap 0.5 " &
    "-a AF_MAX=max\\(AF\\) -a N=first\\(MATCHA_COUNT\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code == 0
  let info = parseInfo(outp)
  # DEL_A_07 has no DB match.
  doAssert "DEL_A_07" in info
  let a07 = info["DEL_A_07"]
  doAssert a07["N"] == "0", "N should be 0 on unmatched, got " & a07["N"]
  doAssert "AF_MAX" notin a07, "AF_MAX should be absent on unmatched"

timed("A19", "end-to-end: INS pass-through; BND annotated via BND_DB"):
  let (outp, code) = run("anno --min-overlap 0.5 -a X=max\\(AF\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code == 0
  let info = parseInfo(outp)
  # INS is still silently dropped at preproc, so the record passes through
  # phase 3 unannotated. BND now flows through matching, so BND_A_11
  # picks up AF from BND_DB_11 (AF=0.40).
  doAssert "INS_A_11" in info, "INS_A_11 should pass through"
  doAssert "X" notin info["INS_A_11"], "INS should not be annotated"
  doAssert "BND_A_11" in info, "BND_A_11 should pass through"
  doAssert "X" in info["BND_A_11"], "BND_A_11 should have AF annotation"
  doAssert info["BND_A_11"]["X"].startsWith("0.4"),
    "BND_A_11 X should be 0.40, got " & info["BND_A_11"]["X"]

timed("A20", "end-to-end: -o file output matches stdout output"):
  let (stdoutText, code1) = run("anno --min-overlap 0.5 -a X=max\\(AF\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code1 == 0
  let tmpVcf = getTempDir() / "matcha_anno_test.vcf"
  defer:
    if fileExists(tmpVcf): removeFile(tmpVcf)
  let (_, code2) = run("anno --min-overlap 0.5 -a X=max\\(AF\\) -o " &
    tmpVcf & " " & FixtureA & " " & FixtureDB)
  doAssert code2 == 0
  doAssert fileExists(tmpVcf)
  doAssert readFile(tmpVcf) == stdoutText

timed("A21", "end-to-end: .bcf output written; --write-index produces .csi"):
  let tmpBcf = getTempDir() / "matcha_anno_test.bcf"
  defer:
    if fileExists(tmpBcf): removeFile(tmpBcf)
    if fileExists(tmpBcf & ".csi"): removeFile(tmpBcf & ".csi")
  let (_, code1) = run("anno --min-overlap 0.5 -a X=max\\(AF\\) -o " &
    tmpBcf & " " & FixtureA & " " & FixtureDB)
  doAssert code1 == 0
  doAssert fileExists(tmpBcf), "bcf not written"
  doAssert not fileExists(tmpBcf & ".csi"),
    "csi should NOT be written by default (--write-index opt-in)"
  let (_, code2) = run("anno --min-overlap 0.5 --write-index -a X=max\\(AF\\) -o " &
    tmpBcf & " " & FixtureA & " " & FixtureDB)
  doAssert code2 == 0
  doAssert fileExists(tmpBcf & ".csi"), "csi not written with --write-index"

timed("A23", "end-to-end: best() ranks by similarity, tie-breaks by posB"):
  # On DEL_A_06 both matches have the same overlap (0.95) → tied; the
  # tiebreak by earliest posB chooses DEL_DB_06b (posB=12950, CALLER=delly).
  let (outp, code) = run("anno --min-overlap 0.5 -a CB=best\\(CALLER\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code == 0
  let info = parseInfo(outp)
  doAssert info["DEL_A_06"]["CB"] == "delly",
    "best(CALLER) on DEL_A_06 should be delly (earliest posB), got " &
    info["DEL_A_06"]["CB"]

timed("A24", "end-to-end: --best-metric is rejected (removed flag)"):
  let (outp, code) = runMerged(
    "anno --min-overlap 0.5 --best-metric overlap -a X=max\\(AF\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code != 0, "--best-metric should be rejected"
  doAssert "best-metric" in outp or "unknown option" in outp,
    "error should mention unknown option, got: " & outp

timed("A25", "end-to-end: BND-DB annotation picks up similarity vector"):
  # BND_A_11 matches BND_DB_11 with sim=1.0.
  let (outp, code) = run("anno --min-overlap 0.5 " &
    "-a S=max\\(MATCHA_SIMILARITY\\) -a N=first\\(MATCHA_COUNT\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code == 0
  let info = parseInfo(outp)
  doAssert "BND_A_11" in info
  doAssert info["BND_A_11"]["N"] == "1",
    "BND_A_11 should have exactly 1 match, got N=" & info["BND_A_11"]["N"]
  let simStr = info["BND_A_11"]["S"]
  doAssert abs(parseFloat(simStr) - 1.0) < 1e-5,
    "BND_A_11 sim should be ~1.0, got " & simStr

timed("A26", "end-to-end: ##matcha_metric header line is emitted"):
  let (outp, code) = run(
    "anno --min-overlap 0.5 -a X=max\\(AF\\) " & FixtureA & " " & FixtureDB)
  doAssert code == 0
  doAssert "##matcha_metric=overlap" in outp,
    "expected ##matcha_metric=overlap in header"
  let (outp2, code2) = run(
    "anno --min-jaccard 0.5 -a X=max\\(AF\\) " & FixtureA & " " & FixtureDB)
  doAssert code2 == 0
  doAssert "##matcha_metric=jaccard" in outp2,
    "expected ##matcha_metric=jaccard in header"

timed("A27", "end-to-end: MATCHA_OVERLAP/MATCHA_JACCARD rejected as unknown DB fields"):
  # The old names are gone — they fall through to the standard "unknown
  # DB SRCFIELD" path.
  let (outp, code) = runMerged(
    "anno --min-overlap 0.5 -a X=max\\(MATCHA_JACCARD\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code != 0
  doAssert "MATCHA_JACCARD" in outp,
    "error should mention the unknown SRCFIELD, got: " & outp

timed("A28", "end-to-end: OUTFIELD collision with input header rejected without --overwrite"):
  # SVTYPE is already in fixtureA's header. Trying to write to it without
  # --overwrite must error.
  let (outp, code) = runMerged("anno --min-overlap 0.5 -a SVTYPE=first\\(CALLER\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code != 0
  doAssert "overwrite" in outp.toLowerAscii or "already" in outp.toLowerAscii,
    "expected overwrite-required error, got: " & outp

timed("A29", "end-to-end: --overwrite allows replacing existing INFO field"):
  let (_, code) = run("anno --overwrite --min-overlap 0.5 -a SVTYPE=first\\(CALLER\\) " &
    FixtureA & " " & FixtureDB)
  doAssert code == 0, "should succeed with --overwrite, exit " & $code
