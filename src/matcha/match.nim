## match.nim — match-mode driver: per-job pair streaming, thread pool,
## main-thread chr:pos resolution, and TSV output.
##
## Workers return MatchPair triples (no resolution). The main thread opens
## slim BCF handles from fileList, queries each representative position via
## CSI `chrom:pos-pos`, and writes TSV rows directly — no MatchResult type.

import std/[atomics, strutils]
import hts
import utils, preproc, log, matchcore

# ---------------------------------------------------------------------------
# Chr:pos slim-BCF resolution (main thread only)
# ---------------------------------------------------------------------------

proc resolveRecord(vcf: VCF; chrom: string; pos, srcIndex: int32;
                   bnd: bool; idxScratch, endScratch: var seq[int32]
                  ): tuple[id: string; endP: int64] =
  ## CSI query `chrom:pos-pos`; return (ID, END) for the record whose
  ## SRC_INDEX matches `srcIndex`. END is 0 for BND (rendered as ".").
  for v in vcf.query(chrom & ":" & $pos & "-" & $pos):
    if readSrcIndex(v, idxScratch) != srcIndex: continue
    var endP: int64
    if not bnd:
      discard extractEnd(v, endScratch, endScratch, endP)
    return (id: $v.ID, endP: endP)

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
      logV("job " & j.chrom & "/" & $j.svtype & "/bin" & $j.binA &
           ": " & $state.results[idx].len & " " & state.label)

proc runJobsWithPool[R](jobs: seq[MatchJob]; cfg: MatchConfig;
                        runner: JobRunner[R]; label: string): seq[seq[R]] =
  result = newSeq[seq[R]](jobs.len)
  if jobs.len == 0: return
  if cfg.nThreads == 1:
    for i, job in jobs.pairs:
      result[i] = runner(job, cfg)
      logV("job " & job.chrom & "/" & $job.svtype & "/bin" & $job.binA &
           ": " & $result[i].len & " " & label)
  else:
    logV("starting " & $cfg.nThreads & " worker thread(s) for " &
         $jobs.len & " job(s)")
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
    logV(label & " workers complete")

proc dispatchPairJob*(job: MatchJob, cfg: MatchConfig): seq[MatchPair] =
  ## Stream MatchPairs for one job. Used by match, collapse, and anno.
  if job.svtype == svBND: streamBndJobPairs(job, cfg)
  else:                   streamJobPairs(job, cfg)

proc runMatchPairJobsWithPool*(jobs: seq[MatchJob]; cfg: MatchConfig): seq[seq[MatchPair]] =
  runJobsWithPool[MatchPair](jobs, cfg, dispatchPairJob, "pairs")

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

proc runMatch*(cfg: MatchConfig) =
  logV("matcha match: A=" & cfg.callsetA & " B=" & cfg.callsetB &
       " threads=" & $cfg.nThreads & " tmp=" & cfg.tmpDir)

  var filesA, filesB: PreprocOutput
  if cfg.selfMode:
    logV("self mode: preprocessing single input")
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A",
                           ioThreads = if cfg.nThreads >= 2: 2 else: 0)
    filesB = filesA
  elif cfg.nThreads >= 2:
    logV("preprocessing A and B in parallel")
    (filesA, filesB) = runParallelPreproc(
      PreprocInput(path: cfg.callsetA, tmpDir: cfg.tmpDir, prefix: "A", ioThreads: 2),
      PreprocInput(path: cfg.callsetB, tmpDir: cfg.tmpDir, prefix: "B", ioThreads: 2))
  else:
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A")
    filesB = preprocessVcf(cfg.callsetB, cfg.tmpDir, "B")

  let (jobs, fileList) = buildWorkQueue(filesA, filesB, cfg)
  logV("work queue: " & $jobs.len & " (chrom, svtype, binA) job(s)")
  if jobs.len == 0:
    return

  let jobResults = runMatchPairJobsWithPool(jobs, cfg)

  # Open slim BCF handles lazily, indexed by position in fileList.
  var slimHandles = newSeq[VCF](fileList.len)
  for i, path in fileList:
    if not open(slimHandles[i], path):
      raise newException(IOError, "cannot open slim BCF: " & path)

  let chromOrder = filesA.chromOrder
  let outFile =
    if isStdoutPath(cfg.outputPath): stdout
    else: open(cfg.outputPath, fmWrite)

  outFile.writeLine("##matcha_metric=" & $cfg.metric)
  outFile.writeLine(OutputHeader)

  var totalMatches = 0
  var idxScratch, endScratch: seq[int32]

  for jrs in jobResults:
    for p in jrs:
      if p.srcIndexB == NO_MATCH: continue   # skip singletons
      let bnd = SvType(p.svtype) == svBND
      let chrom = chromOrder[p.chromIdx]
      let fA = resolveRecord(slimHandles[p.fileIdxA], chrom, p.posA, p.srcIndexA,
                             bnd, idxScratch, endScratch)
      let fB = resolveRecord(slimHandles[p.fileIdxB], chrom, p.posB, p.srcIndexB,
                             bnd, idxScratch, endScratch)
      let endAStr = if bnd: "." else: $fA.endP
      let endBStr = if bnd: "." else: $fB.endP
      outFile.writeLine(
        chrom & "\t" & $p.posA & "\t" & endAStr & "\t" & fA.id & "\t" &
        chrom & "\t" & $p.posB & "\t" & endBStr & "\t" & fB.id & "\t" &
        $SvType(p.svtype) & "\t" & formatFloat(float64(p.sim), ffDecimal, 6))
      inc totalMatches

  for h in slimHandles.mitems: h.close()

  if not isStdoutPath(cfg.outputPath):
    outFile.close()
  logV("wrote " & $totalMatches & " match(es)" &
       (if not isStdoutPath(cfg.outputPath): " to " & cfg.outputPath else: " to stdout"))

  removeTempBcfs(filesA.paths)
  if not cfg.selfMode:
    removeTempBcfs(filesB.paths)
