## match.nim — match-mode driver: per-job pair streaming, thread pool,
## main-thread chr:pos resolution, and TSV output.
##
## Workers return MatchPair triples (no resolution). The main thread opens
## slim BCF handles from fileList, queries each representative position via
## CSI `chrom:pos-pos`, and writes TSV rows directly — no MatchResult type.

import std/[atomics, os, sequtils, sets, strutils]
import hts
import utils, preproc, log, matchcore

# ---------------------------------------------------------------------------
# Chr:pos slim-BCF resolution (main thread only)
# ---------------------------------------------------------------------------

proc formatExtraInfo(v: Variant; fields: seq[string];
                     iBuf: var seq[int32]; fBuf: var seq[float32];
                     sBuf: var string): string =
  ## Emit requested INFO fields as semicolon-delimited KEY=VALUE pairs (VCF style).
  ## Absent fields are silently omitted; returns "." when nothing is present.
  var parts: seq[string]
  for name in fields:
    if v.info().get(name, iBuf) == Status.OK and iBuf.len > 0:
      parts.add(name & "=" & iBuf.mapIt($it).join(","))
    elif v.info().get(name, fBuf) == Status.OK and fBuf.len > 0:
      parts.add(name & "=" & fBuf.mapIt(formatFloat(float64(it), ffDecimal, 6)).join(","))
    elif v.info().get(name, sBuf) == Status.OK and sBuf.len > 0:
      parts.add(name & "=" & sBuf)
  if parts.len == 0: "." else: parts.join(";")

proc resolveRecord(vcf: VCF; chrom: string; pos, srcIndex: int32;
                   infoFields: seq[string]; idxScratch: var seq[int32];
                   iBuf: var seq[int32]; fBuf: var seq[float32]; sBuf: var string
                  ): tuple[id: string; infoStr: string] =
  ## CSI query `chrom:pos-pos`; return (ID, INFO_str) for the record whose
  ## SRC_INDEX matches `srcIndex`.
  for v in vcf.query(chrom & ":" & $pos & "-" & $pos):
    if readSrcIndex(v, idxScratch) != srcIndex: continue
    let infoStr = if infoFields.len > 0: formatExtraInfo(v, infoFields, iBuf, fBuf, sBuf)
                  else: ""
    return (id: $v.ID, infoStr: infoStr)

# ---------------------------------------------------------------------------
# Generic per-job thread pool
# ---------------------------------------------------------------------------

type
  JobRunner[R] = proc(job: MatchJob, cfg: MatchConfig): seq[R] {.nimcall.}

  PoolState[R] = object
    jobs:    seq[MatchJob]
    cfg:     MatchConfig
    next:    Atomic[int]
    results: seq[seq[R]]
    runner:  JobRunner[R]
    label:   string

proc poolWorker[R](state: ptr PoolState[R]) {.thread.} =
  {.cast(gcsafe).}:
    while true:
      let idx = state.next.fetchAdd(1, moRelaxed)
      if idx >= state.jobs.len: break
      state.results[idx] = state.runner(state.jobs[idx], state.cfg)
      let j = state.jobs[idx]
      logVerbose("job " & j.chrom & "/" & $j.svtype & "/bin" & $j.binA &
               ": " & $state.results[idx].len & " " & state.label)

proc runJobsWithPool[R](jobs: seq[MatchJob]; cfg: MatchConfig;
                        runner: JobRunner[R]; label: string;
                        tag = ""): seq[seq[R]] =
  result = newSeq[seq[R]](jobs.len)
  if jobs.len == 0: return
  if cfg.nThreads == 1:
    for i, job in jobs.pairs:
      result[i] = runner(job, cfg)
      logVerbose("job " & job.chrom & "/" & $job.svtype & "/bin" & $job.binA &
               ": " & $result[i].len & " " & label)
  else:
    let suffix = if tag.len > 0: " (" & tag & ")" else: ""
    logInfo("starting " & $cfg.nThreads & " worker thread(s) for " &
            $jobs.len & " job(s)" & suffix)
    var state = PoolState[R](jobs: jobs, cfg: cfg,
                             results: newSeq[seq[R]](jobs.len),
                             runner: runner, label: label)
    state.next.store(0, moRelaxed)
    var threads = newSeq[Thread[ptr PoolState[R]]](cfg.nThreads)
    for i in 0 ..< cfg.nThreads:
      createThread(threads[i], poolWorker[R], addr state)
    for i in 0 ..< cfg.nThreads:
      joinThread(threads[i])
    result = state.results
    logVerbose(label & " workers complete")

proc dispatchPairJob*(job: MatchJob, cfg: MatchConfig): seq[MatchPair] =
  ## Stream MatchPairs for one job. Used by match, collapse, and anno.
  if   job.svtype == svBND: streamBndJobPairs(job, cfg)
  elif job.svtype == svINS: streamInsJobPairs(job, cfg)
  else:                     streamJobPairs(job, cfg)

proc runMatchPairJobsWithPool*(jobs: seq[MatchJob]; cfg: MatchConfig;
                               tag = ""): seq[seq[MatchPair]] =
  runJobsWithPool[MatchPair](jobs, cfg, dispatchPairJob, "pairs", tag)

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

proc runMatch*(cfg: MatchConfig) =
  logInfo("matcha match: A=" & cfg.callsetA & " B=" & cfg.callsetB &
          " threads=" & $cfg.nThreads & " tmp=" & cfg.tmpDir)

  let extra = cfg.infoFields   # shorthand
  let keptChrsSet = toHashSet(cfg.keptChrs)
  let chrSetSet   = toHashSet(cfg.chrSet)

  var filesA, filesB: PreprocOutput
  if cfg.selfMode:
    logInfo("self mode: preprocessing single input")
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A", extra,
                           ioThreads = if cfg.nThreads >= 2: 2 else: 0,
                           keptChrs = keptChrsSet, chrSet = chrSetSet)
    filesB = filesA
  elif cfg.nThreads >= 2:
    logInfo("preprocessing A and B in parallel")
    (filesA, filesB) = runParallelPreproc(
      PreprocInput(path: cfg.callsetA, tmpDir: cfg.tmpDir, prefix: "A",
                   extraKeep: extra, ioThreads: 2,
                   keptChrs: keptChrsSet, chrSet: chrSetSet),
      PreprocInput(path: cfg.callsetB, tmpDir: cfg.tmpDir, prefix: "B",
                   extraKeep: extra, ioThreads: 2,
                   keptChrs: keptChrsSet, chrSet: chrSetSet))
  else:
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A", extra,
                           keptChrs = keptChrsSet, chrSet = chrSetSet)
    filesB = preprocessVcf(cfg.callsetB, cfg.tmpDir, "B", extra,
                           keptChrs = keptChrsSet, chrSet = chrSetSet)
  var chrsSeen = filesA.chrsSeen
  for c in filesB.chrsSeen: chrsSeen.incl(c)
  warnMissingChrs(cfg.keptChrs, chrsSeen)

  let (jobs, fileList) = buildWorkQueue(filesA, filesB, cfg)
  logInfo("work queue: " & $jobs.len & " (chrom, svtype, binA) job(s)")
  if jobs.len == 0:
    return

  let jobResults = runMatchPairJobsWithPool(jobs, cfg)

  var slimHandles = newSeq[VCF](fileList.len)
  for i, path in fileList:
    if not open(slimHandles[i], path):
      raise newException(IOError, "cannot open slim BCF: " & path)

  let chromOrder = filesA.chromOrder
  let outFile =
    if isStdoutPath(cfg.outputPath): stdout
    else: open(cfg.outputPath, fmWrite)

  outFile.writeLine("##matcha_metric=" & $cfg.metric)
  if extra.len == 0:
    outFile.writeLine(OutputHeader)
  else:
    outFile.writeLine(
      "#CHROM_A\tPOS_A\tID_A\tINFO_A\tCHROM_B\tPOS_B\tID_B\tINFO_B\tSVTYPE\tSIMILARITY")

  var totalMatches = 0
  var idxScratch, iBuf: seq[int32]
  var fBuf: seq[float32]
  var sBuf: string

  for jrs in jobResults:
    for p in jrs:
      if p.srcIndexB == NO_MATCH: continue
      let chrom  = chromOrder[p.chromIdx]
      let simStr = formatFloat(float64(p.sim), ffDecimal, 6)
      let svStr  = $SvType(p.svtype)
      let fA = resolveRecord(slimHandles[p.fileIdxA], chrom, p.posA, p.srcIndexA,
                             extra, idxScratch, iBuf, fBuf, sBuf)
      let fB = resolveRecord(slimHandles[p.fileIdxB], chrom, p.posB, p.srcIndexB,
                             extra, idxScratch, iBuf, fBuf, sBuf)
      if extra.len == 0:
        outFile.writeLine(
          chrom & "\t" & $p.posA & "\t" & fA.id & "\t" &
          chrom & "\t" & $p.posB & "\t" & fB.id & "\t" &
          svStr & "\t" & simStr)
      else:
        outFile.writeLine(
          chrom & "\t" & $p.posA & "\t" & fA.id & "\t" & fA.infoStr & "\t" &
          chrom & "\t" & $p.posB & "\t" & fB.id & "\t" & fB.infoStr & "\t" &
          svStr & "\t" & simStr)
      inc totalMatches

  for h in slimHandles.mitems: h.close()

  if not isStdoutPath(cfg.outputPath):
    outFile.close()
  logInfo("wrote " & $totalMatches & " match(es)" &
          (if not isStdoutPath(cfg.outputPath): " to " & cfg.outputPath else: " to stdout"))

  removeDir(cfg.tmpDir)
