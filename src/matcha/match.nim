## match.nim — match-mode driver: per-job matching adapter over matchcore,
## thread pool, and output assembly.

import std/[atomics, os, tables]
import hts
import utils, preproc, log, bins, matchcore

# ---------------------------------------------------------------------------
# Per-job matching (adapter over the shared streamJobPairs)
# ---------------------------------------------------------------------------

proc runMatchJob*(job: MatchJob, cfg: MatchConfig): seq[MatchResult] =
  ## Stream A records, query each adjacent populated B bin via TiledBuffers,
  ## apply the metric filter, and emit MatchResults. Self-mode dedup happens
  ## in the emit callback: aOff < bOff selects one copy of each symmetric
  ## pair and excludes the trivial X-vs-X case.
  let selfMode = cfg.selfMode
  let svtype = job.svtype
  streamJobPairs[bool, MatchResult](job, cfg,
    extract = proc(v: Variant): bool = false,
    emit = proc(va: Variant; posA, endA, aOff: int64;
                cand: BufferedRec; bExtra: bool;
                ovl, jac: float64): PairResult[MatchResult] =
      if selfMode and aOff >= cand.bOffset:
        return PairResult[MatchResult](keep: false)
      PairResult[MatchResult](keep: true, item: MatchResult(
        chrom:   $va.CHROM,
        posA:    posA, endA: endA, idA: $va.ID,
        posB:    cand.pos, endB: cand.endPos, idB: cand.id,
        svtype:  svtype,
        overlap: ovl, jaccard: jac,
        aOffset: aOff, bOffset: cand.bOffset,
      )))

# ---------------------------------------------------------------------------
# Thread pool (global state + atomic counter + per-slot results)
# ---------------------------------------------------------------------------

var gJobs:        seq[MatchJob]
var gCfg:         MatchConfig
var gNextJob:     Atomic[int]
var gJobResults:  seq[seq[MatchResult]]   # disjoint slot per job index

proc threadWorker(dummy: int) {.thread.} =
  {.cast(gcsafe).}:
    while true:
      let idx = gNextJob.fetchAdd(1, moRelaxed)
      if idx >= gJobs.len:
        break
      gJobResults[idx] = runMatchJob(gJobs[idx], gCfg)
      logV("job " & gJobs[idx].chrom & "/" & $gJobs[idx].svtype &
           "/bin" & $gJobs[idx].binA &
           ": " & $gJobResults[idx].len & " matches")

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
    filesB = filesA   # same paths, bins, chroms — dedup happens later
  elif cfg.nThreads >= 2:
    logV("preprocessing A and B in parallel")
    (filesA, filesB) = runParallelPreproc(
      PreprocInput(path: cfg.callsetA, tmpDir: cfg.tmpDir, prefix: "A", ioThreads: 2),
      PreprocInput(path: cfg.callsetB, tmpDir: cfg.tmpDir, prefix: "B", ioThreads: 2))
  else:
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A")
    filesB = preprocessVcf(cfg.callsetB, cfg.tmpDir, "B")

  var jobs = buildWorkQueue(filesA, filesB, cfg)
  logV("work queue: " & $jobs.len & " (chrom, svtype, binA) job(s)")
  if jobs.len == 0:
    return

  var jobResults = newSeq[seq[MatchResult]](jobs.len)
  if cfg.nThreads == 1:
    for i, job in jobs.pairs:
      jobResults[i] = runMatchJob(job, cfg)
      logV("job " & job.chrom & "/" & $job.svtype & "/bin" & $job.binA &
           ": " & $jobResults[i].len & " matches")
  else:
    logV("starting " & $cfg.nThreads & " worker thread(s)")
    gJobs = jobs
    gCfg  = cfg
    gJobResults = newSeq[seq[MatchResult]](jobs.len)
    gNextJob.store(0, moRelaxed)
    var threads = newSeq[Thread[int]](cfg.nThreads)
    for i in 0 ..< cfg.nThreads:
      createThread(threads[i], threadWorker, i)
    for i in 0 ..< cfg.nThreads:
      joinThread(threads[i])
    jobResults = gJobResults
    logV("workers complete")

  # Write the header line plus all matches in deterministic (job-sorted) order.
  let outFile =
    if cfg.outputPath == "": stdout
    else: open(cfg.outputPath, fmWrite)

  outFile.writeLine(OutputHeader)

  var totalMatches = 0
  for jrs in jobResults:
    for mr in jrs:
      outFile.writeLine(formatMatchResult(mr))
      inc totalMatches

  if cfg.outputPath != "":
    outFile.close()
  logV("wrote " & $totalMatches & " match(es)" &
       (if cfg.outputPath != "": " to " & cfg.outputPath else: " to stdout"))

  # Clean up temp BCFs (per-(svtype, bin) files for both callsets).
  # In self mode, filesB aliases filesA — skip the second loop.
  for path in filesA.paths.values:
    if fileExists(path):     removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")
  if not cfg.selfMode:
    for path in filesB.paths.values:
      if fileExists(path):     removeFile(path)
      if fileExists(path & ".csi"): removeFile(path & ".csi")
