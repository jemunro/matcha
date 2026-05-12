## anno.nim — `matcha anno` subcommand.
##
## Annotate input.vcf with INFO fields aggregated from database.vcf, based on
## SV matches computed by the shared binning + tiled-buffer engine.
##
## Three phases:
##   1. Preproc A (default keep-set) and B (default keep-set ∪ user SRCFIELDs).
##   2. Match A vs B via runAnnoJob — extends runMatchJob to also pull the
##      user-requested INFO values off each B candidate during the same pass.
##   3. Stream the *original* A file, capture each record's bgzf_tell offset,
##      look up its match set in an in-memory map keyed by aOffset, apply
##      aggregation functions, and write the annotated record.

import std/[algorithm, atomics, os, sequtils, sets, strutils, tables]
import hts
import hts/private/hts_concat
import utils, preproc, log, bins, matchcore

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

type
  AggFunc* = enum
    afMax     = "max"
    afMin     = "min"
    afMean    = "mean"
    afFirst   = "first"
    afLast    = "last"
    afBest    = "best"
    afAll     = "all"
    afUnique  = "unique"

  MatchaVar* = enum
    ## Which (if any) implicit MATCHA variable the SRCFIELD refers to.
    ## MATCHA_SIMILARITY collapses what used to be MATCHA_OVERLAP /
    ## MATCHA_JACCARD into a single per-match Float vector — the active
    ## metric for interval rows, the slop-based proximity for BND rows.
    mvNone, mvCount, mvSimilarity

  AnnoExpr* = object
    outField*:   string        ## INFO field name to emit
    fn*:         AggFunc
    srcField*:   string        ## raw token after FUNC( ... )
    matchaVar*:  MatchaVar     ## mvNone unless srcField names a MATCHA_*
    dbType*:     string        ## "Integer" / "Float" / "String" (DB-sourced only)
    dbNumber*:   string        ## "1" / "." / "N" (DB-sourced only)

  AnnoConfig* = object
    metric*:        Metric       ## Active interval metric (mOverlap | mJaccard).
    threshold*:     float64      ## Minimum score for the active metric.
    bndSlop*:       int          ## --bnd-slop (default 100)
    overwrite*:     bool
    nThreads*:      int
    tmpDir*:        string
    outputPath*:    string
    callsetA*:      string     ## input
    callsetB*:      string     ## database
    exprs*:         seq[AnnoExpr]
    dbFields*:      seq[string]  ## deduped DB SRCFIELDs (excluding MATCHA_*)

  AnnoMatch* = object
    aOffset*:    int64
    bOffset*:    int64
    posB*:       int64
    similarity*: float64
    payload*:    Table[string, seq[string]]  ## DB srcField → stringified values

  DbPayload* = Table[string, seq[string]]

# ---------------------------------------------------------------------------
# Expression parsing and validation
# ---------------------------------------------------------------------------

const ValidIdChars = {'A'..'Z', 'a'..'z', '0'..'9', '_'}

proc isValidInfoName(s: string): bool =
  if s.len == 0: return false
  for c in s:
    if c notin ValidIdChars: return false
  true

proc parseAnnoExpr*(s: string): AnnoExpr =
  ## Parse `OUTFIELD=FUNC(SRCFIELD)`. Raises ValueError on malformed input.
  let eqIdx = s.find('=')
  if eqIdx < 1:
    raise newException(ValueError,
      "invalid -a expression (need OUTFIELD=FUNC(SRCFIELD)): " & s)
  result.outField = s[0 ..< eqIdx].strip
  if not isValidInfoName(result.outField):
    raise newException(ValueError,
      "invalid OUTFIELD name '" & result.outField & "' in expression: " & s)

  let rest = s[eqIdx + 1 .. ^1].strip
  let lp = rest.find('(')
  if lp < 1 or not rest.endsWith(")"):
    raise newException(ValueError,
      "invalid -a expression (need FUNC(SRCFIELD)): " & s)
  let fnStr = rest[0 ..< lp].strip.toLowerAscii
  let srcStr = rest[lp + 1 .. ^2].strip
  result.fn =
    case fnStr
    of "max":    afMax
    of "min":    afMin
    of "mean":   afMean
    of "first":  afFirst
    of "last":   afLast
    of "best":   afBest
    of "all":    afAll
    of "unique": afUnique
    else:
      raise newException(ValueError,
        "unknown aggregation function '" & fnStr & "' in expression: " & s)
  if not isValidInfoName(srcStr):
    raise newException(ValueError,
      "invalid SRCFIELD name '" & srcStr & "' in expression: " & s)
  result.srcField = srcStr
  result.matchaVar =
    case srcStr
    of "MATCHA_COUNT":      mvCount
    of "MATCHA_SIMILARITY": mvSimilarity
    else:                   mvNone

proc outputNumberType*(e: AnnoExpr): tuple[number, typ: string] =
  ## Derive the output INFO Number/Type from the expression's FUNC and source.
  case e.matchaVar
  of mvCount:
    return ("1", "Integer")
  of mvSimilarity:
    case e.fn
    of afAll, afUnique: return (".", "Float")
    else:               return ("1", "Float")
  of mvNone:
    case e.fn
    of afMean:           return ("1", "Float")
    of afAll, afUnique:  return (".", e.dbType)
    else:                return ("1", e.dbType)

proc validateAnnoExprs*(cfg: var AnnoConfig) =
  ## Validate `-a` expressions against the DB and input headers. Hard-errors
  ## on duplicate OUTFIELD, unknown DB SRCFIELD, or OUTFIELD collisions in the
  ## input header (when --overwrite not set). Populates dbType/dbNumber on
  ## each non-MATCHA expression and the dbFields list on the config.
  if cfg.exprs.len == 0:
    stderr.writeLine "error: anno requires at least one -a expression"
    quit(1)

  var seenOut = initHashSet[string]()
  for e in cfg.exprs:
    if e.outField in seenOut:
      stderr.writeLine "error: duplicate OUTFIELD '" & e.outField & "' in -a expressions"
      quit(1)
    seenOut.incl(e.outField)

  # DB header lookup for non-MATCHA SRCFIELDs.
  var vcfDB: VCF
  if not open(vcfDB, cfg.callsetB):
    stderr.writeLine "error: cannot open database: " & cfg.callsetB
    quit(1)
  var dbFieldsSet = initHashSet[string]()
  for i in 0 ..< cfg.exprs.len:
    if cfg.exprs[i].matchaVar != mvNone: continue
    let src = cfg.exprs[i].srcField
    try:
      let hr = vcfDB.header.get(src, BCF_HEADER_TYPE.BCF_HL_INFO)
      cfg.exprs[i].dbType = hr["Type"]
      cfg.exprs[i].dbNumber = hr["Number"]
    except KeyError:
      stderr.writeLine "error: SRCFIELD '" & src & "' not found in database header: " & cfg.callsetB
      quit(1)
    dbFieldsSet.incl(src)
  vcfDB.close()
  cfg.dbFields = toSeq(dbFieldsSet.items)
  cfg.dbFields.sort()

  # Input-header OUTFIELD collision check.
  var vcfA: VCF
  if not open(vcfA, cfg.callsetA):
    stderr.writeLine "error: cannot open input: " & cfg.callsetA
    quit(1)
  for e in cfg.exprs:
    try:
      discard vcfA.header.get(e.outField, BCF_HEADER_TYPE.BCF_HL_INFO)
      # already present
      if not cfg.overwrite:
        stderr.writeLine "error: OUTFIELD '" & e.outField &
          "' already in input header; pass --overwrite to replace"
        quit(1)
      else:
        logWarn("OUTFIELD '" & e.outField & "' present in input header — replacing (--overwrite)")
    except KeyError:
      discard
  vcfA.close()

# ---------------------------------------------------------------------------
# B-INFO extraction during matching
# ---------------------------------------------------------------------------

proc extractDbValues(v: Variant; dbFields: seq[string];
                     dbTypes: Table[string, string];
                     iBuf: var seq[int32];
                     fBuf: var seq[float32];
                     sBuf: var string): DbPayload =
  ## Pull each requested DB INFO field off the given variant, stringifying
  ## scalar/list values into a parallel seq[string]. Absent fields are
  ## omitted from the payload; aggregation treats missing entries as
  ## "no contribution from this match".
  for name in dbFields:
    let typ = dbTypes.getOrDefault(name, "String")
    case typ
    of "Integer":
      if v.info().get(name, iBuf) == Status.OK:
        var vals: seq[string]
        for x in iBuf: vals.add($x)
        result[name] = vals
    of "Float":
      if v.info().get(name, fBuf) == Status.OK:
        var vals: seq[string]
        for x in fBuf: vals.add(formatFloat(float64(x), ffDecimal, 6))
        result[name] = vals
    else:
      # String / Flag — htslib returns the raw comma-joined storage.
      if v.info().get(name, sBuf) == Status.OK and sBuf.len > 0:
        result[name] = sBuf.split(',')

proc runAnnoJob*(job: MatchJob; cfg: AnnoConfig;
                 dbTypes: Table[string, string]): seq[AnnoMatch] =
  ## Anno-mode adapter over `streamJobPairs` (interval matches). The
  ## `extract` callback pulls user-requested DB INFO values off each B
  ## candidate during the same tile fetch that match-mode uses; the `emit`
  ## callback wraps the pair into an AnnoMatch. No self-mode filter — anno
  ## is asymmetric (input vs database) by construction.
  let dbFields = cfg.dbFields
  # Scratch buffers reused across all extract() calls within this job.
  var iBuf: seq[int32]
  var fBuf: seq[float32]
  var sBuf: string

  streamJobPairs[DbPayload, AnnoMatch](job, MatchConfig(
      metric: cfg.metric, threshold: cfg.threshold, bndSlop: cfg.bndSlop),
    extract = proc(v: Variant): DbPayload =
      extractDbValues(v, dbFields, dbTypes, iBuf, fBuf, sBuf),
    emit = proc(va: Variant; posA, endA, aOff: int64;
                cand: BufferedRec; bExtra: DbPayload;
                sim: float64): PairResult[AnnoMatch] =
      PairResult[AnnoMatch](keep: true, item: AnnoMatch(
        aOffset:    aOff, bOffset: cand.bOffset, posB: cand.pos,
        similarity: sim,
        payload:    bExtra,
      )))

proc runAnnoBndJob*(job: MatchJob; cfg: AnnoConfig;
                    dbTypes: Table[string, string]): seq[AnnoMatch] =
  ## Anno-mode BND adapter. Uses streamBndJobPairs with slop-based proximity.
  let dbFields = cfg.dbFields
  var iBuf: seq[int32]
  var fBuf: seq[float32]
  var sBuf: string

  streamBndJobPairs[DbPayload, AnnoMatch](job, MatchConfig(
      metric: cfg.metric, threshold: cfg.threshold, bndSlop: cfg.bndSlop),
    extract = proc(v: Variant): DbPayload =
      extractDbValues(v, dbFields, dbTypes, iBuf, fBuf, sBuf),
    emit = proc(va: Variant; posA, pos2A, aOff: int64;
                cand: BufferedRec; bExtra: DbPayload;
                sim: float64): PairResult[AnnoMatch] =
      PairResult[AnnoMatch](keep: true, item: AnnoMatch(
        aOffset:    aOff, bOffset: cand.bOffset, posB: cand.pos,
        similarity: sim,
        payload:    bExtra,
      )))

proc dispatchAnnoJob(job: MatchJob; cfg: AnnoConfig;
                     dbTypes: Table[string, string]): seq[AnnoMatch] {.inline.} =
  if job.svtype == svBND: runAnnoBndJob(job, cfg, dbTypes)
  else:                   runAnnoJob(job, cfg, dbTypes)

# ---------------------------------------------------------------------------
# Aggregation kernel
# ---------------------------------------------------------------------------

proc parseFloatSafe(s: string): float64 =
  try: parseFloat(s) except ValueError: 0.0

proc selectBestIdx(matches: seq[AnnoMatch]): int =
  ## Index of the match with the highest similarity; ties broken by smaller
  ## posB (earliest by position).
  result = 0
  for i in 1 ..< matches.len:
    if matches[i].similarity > matches[result].similarity:
      result = i
    elif matches[i].similarity == matches[result].similarity and
         matches[i].posB < matches[result].posB:
      result = i

proc gatherValues(e: AnnoExpr; matches: seq[AnnoMatch]): seq[string] =
  ## Build the pooled per-match value list for the expression's SRCFIELD.
  ## For MATCHA_* this is computed from the match metrics; for DB fields
  ## it's pulled from each match's payload (already split on ',' upstream).
  case e.matchaVar
  of mvCount:
    result.add($matches.len)
  of mvSimilarity:
    for m in matches: result.add(formatFloat(m.similarity, ffDecimal, 6))
  of mvNone:
    for m in matches:
      if e.srcField in m.payload:
        for v in m.payload[e.srcField]:
          result.add(v)

proc applyAggFunc*(e: AnnoExpr; matches: seq[AnnoMatch]): seq[string] =
  ## Run the aggregation function over the match set. Returns the formatted
  ## value(s) to set on the output INFO field. Empty seq → field absent.
  if matches.len == 0:
    # MATCHA_COUNT degenerates to 0; everything else is absent.
    if e.matchaVar == mvCount: return @["0"]
    return @[]

  # Sort by posB for deterministic first/last/best tiebreaks.
  var sorted = matches
  sorted.sort(proc(a, b: AnnoMatch): int = cmp(a.posB, b.posB))

  let pooled = gatherValues(e, sorted)

  case e.fn
  of afMax:
    if pooled.len == 0: return @[]
    var m = parseFloatSafe(pooled[0])
    for s in pooled[1 .. ^1]:
      let v = parseFloatSafe(s)
      if v > m: m = v
    if e.matchaVar == mvCount or e.dbType == "Integer":
      return @[$int(m)]
    return @[formatFloat(m, ffDecimal, 6)]
  of afMin:
    if pooled.len == 0: return @[]
    var m = parseFloatSafe(pooled[0])
    for s in pooled[1 .. ^1]:
      let v = parseFloatSafe(s)
      if v < m: m = v
    if e.matchaVar == mvCount or e.dbType == "Integer":
      return @[$int(m)]
    return @[formatFloat(m, ffDecimal, 6)]
  of afMean:
    if pooled.len == 0: return @[]
    var s = 0.0
    for v in pooled: s += parseFloatSafe(v)
    return @[formatFloat(s / float(pooled.len), ffDecimal, 6)]
  of afFirst:
    if pooled.len == 0: return @[]
    return @[pooled[0]]
  of afLast:
    if pooled.len == 0: return @[]
    return @[pooled[^1]]
  of afBest:
    let bi = selectBestIdx(sorted)
    case e.matchaVar
    of mvCount:      return @[$sorted.len]
    of mvSimilarity: return @[formatFloat(sorted[bi].similarity, ffDecimal, 6)]
    of mvNone:
      if e.srcField in sorted[bi].payload and
         sorted[bi].payload[e.srcField].len > 0:
        return @[sorted[bi].payload[e.srcField][0]]
      return @[]
  of afAll:
    return pooled
  of afUnique:
    var seen = initHashSet[string]()
    var uniq: seq[string]
    for v in pooled:
      if v notin seen:
        seen.incl(v); uniq.add(v)
    return uniq

# ---------------------------------------------------------------------------
# Thread pool for matching phase
# ---------------------------------------------------------------------------

var gAnnoJobs:     seq[MatchJob]
var gAnnoCfg:      AnnoConfig
var gAnnoDbTypes:  Table[string, string]
var gAnnoNext:     Atomic[int]
var gAnnoResults:  seq[seq[AnnoMatch]]

proc annoThreadWorker(dummy: int) {.thread.} =
  {.cast(gcsafe).}:
    while true:
      let idx = gAnnoNext.fetchAdd(1, moRelaxed)
      if idx >= gAnnoJobs.len: break
      gAnnoResults[idx] = dispatchAnnoJob(gAnnoJobs[idx], gAnnoCfg, gAnnoDbTypes)
      logV("anno job " & gAnnoJobs[idx].chrom & "/" & $gAnnoJobs[idx].svtype &
           "/bin" & $gAnnoJobs[idx].binA &
           ": " & $gAnnoResults[idx].len & " match(es)")

# ---------------------------------------------------------------------------
# Output assembly (phase 3)
# ---------------------------------------------------------------------------

proc isStdoutPath(p: string): bool = p == "" or p == "-"

proc validateOutputPath(p: string) =
  if isStdoutPath(p): return
  if not (p.endsWith(".vcf") or p.endsWith(".vcf.gz") or p.endsWith(".bcf")):
    stderr.writeLine "error: -o must end in .vcf, .vcf.gz, or .bcf: " & p
    quit(1)

proc registerOutputHeader(h: vcf.Header; cfg: AnnoConfig) =
  # Record the active interval metric. BND rows always use slop-similarity;
  # this line documents which metric drove interval matching in this run.
  let metricName = $cfg.metric
  discard h.add_string("##matcha_metric=" & metricName)
  for e in cfg.exprs:
    let (num, typ) = outputNumberType(e)
    let desc = "matcha anno: " & $e.fn & "(" & e.srcField & ")"
    # If --overwrite and the field already exists, remove it first so the
    # new declaration sticks (htslib's bcf_hdr_append silently no-ops on
    # duplicates).
    try:
      discard h.get(e.outField, BCF_HEADER_TYPE.BCF_HL_INFO)
      if cfg.overwrite:
        discard h.remove_info(e.outField)
    except KeyError:
      discard
    discard h.add_info(e.outField, num, typ, desc)

proc setInfoValues(v: Variant; e: AnnoExpr; values: seq[string]) =
  ## Write one expression's aggregated values onto the variant. Skips when
  ## values is empty (the field stays absent).
  if values.len == 0: return
  let (num, typ) = outputNumberType(e)
  case typ
  of "Integer":
    var ints: seq[int32]
    for s in values:
      try: ints.add(int32(parseInt(s)))
      except ValueError: ints.add(int32(parseFloatSafe(s)))
    if ints.len > 0:
      if num == "1":
        var x = ints[0]
        discard v.info.set(e.outField, x)
      else:
        discard v.info.set(e.outField, ints)
  of "Float":
    var floats: seq[float32]
    for s in values: floats.add(float32(parseFloatSafe(s)))
    if floats.len > 0:
      if num == "1":
        var x = floats[0]
        discard v.info.set(e.outField, x)
      else:
        discard v.info.set(e.outField, floats)
  else:
    # String. htslib stores Number=. String as a single comma-joined string.
    var s =
      if num == "1": values[0]
      else: values.join(",")
    discard v.info.set(e.outField, s)

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

proc runAnno*(cfg: var AnnoConfig) =
  logV("matcha anno: A=" & cfg.callsetA & " B=" & cfg.callsetB &
       " threads=" & $cfg.nThreads & " tmp=" & cfg.tmpDir)
  validateOutputPath(cfg.outputPath)
  validateAnnoExprs(cfg)

  # --- Phase 1: preproc -----------------------------------------------------
  var filesA, filesB: PreprocOutput
  if cfg.nThreads >= 2:
    logV("preprocessing A and B in parallel" &
         (if cfg.dbFields.len > 0: " (B extra: " & cfg.dbFields.join(",") & ")" else: ""))
    (filesA, filesB) = runParallelPreproc(
      PreprocInput(path: cfg.callsetA, tmpDir: cfg.tmpDir, prefix: "A", ioThreads: 2),
      PreprocInput(path: cfg.callsetB, tmpDir: cfg.tmpDir, prefix: "B",
                   extraKeep: cfg.dbFields, ioThreads: 2))
  else:
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A")
    filesB = preprocessVcf(cfg.callsetB, cfg.tmpDir, "B", cfg.dbFields)

  # buildWorkQueue takes a MatchConfig; build a shim with matching thresholds.
  var matchCfg = MatchConfig(
    metric: cfg.metric, threshold: cfg.threshold,
    nThreads: cfg.nThreads, tmpDir: cfg.tmpDir, selfMode: false,
  )
  let jobs = buildWorkQueue(filesA, filesB, matchCfg)
  logV("anno work queue: " & $jobs.len & " (chrom, svtype, binA) job(s)")

  # --- Phase 2: matching with B-INFO extraction -----------------------------
  var dbTypes = initTable[string, string]()
  for e in cfg.exprs:
    if e.matchaVar == mvNone:
      dbTypes[e.srcField] = e.dbType

  var jobResults = newSeq[seq[AnnoMatch]](jobs.len)
  if jobs.len > 0:
    if cfg.nThreads == 1:
      for i, j in jobs.pairs:
        jobResults[i] = dispatchAnnoJob(j, cfg, dbTypes)
        logV("anno job " & j.chrom & "/" & $j.svtype & "/bin" & $j.binA &
             ": " & $jobResults[i].len & " match(es)")
    else:
      gAnnoJobs = jobs
      gAnnoCfg = cfg
      gAnnoDbTypes = dbTypes
      gAnnoResults = newSeq[seq[AnnoMatch]](jobs.len)
      gAnnoNext.store(0, moRelaxed)
      var threads = newSeq[Thread[int]](cfg.nThreads)
      for i in 0 ..< cfg.nThreads:
        createThread(threads[i], annoThreadWorker, i)
      for i in 0 ..< cfg.nThreads:
        joinThread(threads[i])
      jobResults = gAnnoResults

  # Group matches by aOffset.
  var byAoff = initTable[int64, seq[AnnoMatch]]()
  var totalMatches = 0
  for jrs in jobResults:
    for m in jrs:
      byAoff.mgetOrPut(m.aOffset, @[]).add(m)
      inc totalMatches
  logV("anno: collected " & $totalMatches & " match(es) over " &
       $byAoff.len & " annotated A record(s)")

  # --- Phase 3: stream original A, write annotated output -------------------
  var vcfA: VCF
  if not open(vcfA, cfg.callsetA, threads = if cfg.nThreads >= 2: 2 else: 0):
    raise newException(IOError, "cannot reopen input: " & cfg.callsetA)
  vcfA.set_samples(@["^"])
  # Register new INFO fields on the SOURCE header. info.set() resolves field
  # IDs via the variant's source-header pointer; the writer then copies the
  # already-augmented header so reader and writer share the same ID space.
  registerOutputHeader(vcfA.header, cfg)

  var outWriter: VCF
  # hts-nim's `close()` skips hts_close() for fname=="-" — the BGZF buffer
  # never flushes and the user sees no output. "/dev/stdout" routes through
  # the regular file path, so hts_close → bgzf_close flushes correctly.
  let outPath = if isStdoutPath(cfg.outputPath): "/dev/stdout" else: cfg.outputPath
  if not open(outWriter, outPath, mode = "w"):
    raise newException(IOError, "cannot open output: " & outPath)
  outWriter.copy_header(vcfA.header)
  if not outWriter.write_header():
    raise newException(IOError, "failed to write output header")

  # Enable bgzf compression threads and real-time CSI indexing for compressed
  # outputs. Plain VCF (stdout or .vcf) is neither compressed nor indexed.
  let isBgzf = not isStdoutPath(cfg.outputPath) and
               (outPath.endsWith(".bcf") or outPath.endsWith(".vcf.gz"))
  if cfg.nThreads >= 2 and isBgzf:
    discard bgzf_mt(bgzfHandle(outWriter), 2, 128)
  var outIdx: ptr hts_idx_t = nil
  if isBgzf:
    let headerOff = uint64(bgzf_tell(bgzfHandle(outWriter)))
    outIdx = hts_idx_init(0, HTS_FMT_CSI.cint, headerOff, 14, 5)
    if outIdx == nil:
      raise newException(IOError, "cannot create CSI index for: " & outPath)

  # bgzf_tell prime + update mirrors preproc's offset capture.
  var nextOffset = int64(bgzf_tell(bgzfHandle(vcfA)))
  var nInput = 0
  var nAnnotated = 0
  var endBuf: seq[int32]
  for v in vcfA:
    let recOff = nextOffset
    inc nInput
    if recOff in byAoff:
      let matches = byAoff[recOff]
      for e in cfg.exprs:
        let vals = applyAggFunc(e, matches)
        setInfoValues(v, e, vals)
      inc nAnnotated
    else:
      # No matches: still emit MATCHA_COUNT=0 for any expression that
      # wraps MATCHA_COUNT; other expressions stay absent.
      for e in cfg.exprs:
        if e.matchaVar == mvCount:
          var zero: int32 = 0
          discard v.info.set(e.outField, zero)
    let woff = uint64(bgzf_tell(bgzfHandle(outWriter)))
    if not outWriter.write_variant(v):
      raise newException(IOError, "failed to write variant at " &
        $v.CHROM & ":" & $v.POS)
    if outIdx != nil:
      let endPos =
        if v.info.get("END", endBuf) == Status.OK and endBuf.len > 0:
          int64(endBuf[0])
        else:
          int64(v.c.pos + v.c.rlen)
      discard hts_idx_push(outIdx, v.c.rid, v.c.pos, endPos, woff, 1)
    nextOffset = int64(bgzf_tell(bgzfHandle(vcfA)))

  vcfA.close()
  if outIdx != nil:
    let finalOff = uint64(bgzf_tell(bgzfHandle(outWriter)))
    hts_idx_finish(outIdx, finalOff)
  outWriter.close()
  logV("anno: wrote " & $nInput & " record(s) (" &
       $nAnnotated & " annotated) to " &
       (if isStdoutPath(cfg.outputPath): "stdout" else: cfg.outputPath))
  if outIdx != nil:
    hts_idx_save(outIdx, cfg.outputPath.cstring, HTS_FMT_CSI.cint)
    hts_idx_destroy(outIdx)
    logV("indexed " & cfg.outputPath)

  # Clean up temp BCFs (A and B preproc artifacts).
  for path in filesA.paths.values:
    if fileExists(path):           removeFile(path)
    if fileExists(path & ".csi"):  removeFile(path & ".csi")
  for path in filesB.paths.values:
    if fileExists(path):           removeFile(path)
    if fileExists(path & ".csi"):  removeFile(path & ".csi")
