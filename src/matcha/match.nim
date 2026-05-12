## match.nim — match-mode driver: per-job matching adapter over matchcore,
## per-job slim-BCF resolution into MatchResult, thread pool, output assembly.
##
## matchcore.streamJobPairs / streamBndJobPairs return only (aOff, bOff, sim).
## This adapter resolves the (chrom, posA, endA, idA, posB, endB, idB) fields
## by re-scanning the per-job slim BCFs once each for `job.chrom`, picking
## records up by INFO/MATCHA_BOFF. Slim BCFs are keep-set only (no FORMAT or
## samples) so this pass is cheap. svtype and chrom come from `job` directly.

import std/[atomics, os, sets, tables]
import hts
import utils, preproc, log, matchcore

# ---------------------------------------------------------------------------
# Per-job resolution: slim-BCF MATCHA_BOFF → (POS, END, ID) for output rows
# ---------------------------------------------------------------------------

type SlimFields = tuple[pos, endP: int64; id: string]

proc collectFieldsA(path, chrom: string;
                    needed: HashSet[int64]): Table[int64, SlimFields] =
  ## Scan a slim BCF for `chrom`, capturing (POS, END, ID) keyed by
  ## MATCHA_BOFF for every record whose offset is in `needed`. END is
  ## decoded via the shared extractEnd helper.
  var vcf: VCF
  if not open(vcf, path):
    raise newException(IOError, "cannot reopen slim BCF: " & path)
  var endData, svlenData, boffData: seq[int32]
  for v in vcf.query(chrom):
    let off = readBoff(v, boffData)
    if off notin needed: continue
    var endP: int64
    if not extractEnd(v, endData, svlenData, endP): continue
    result[off] = (v.POS, endP, $v.ID)
  vcf.close()

proc collectFieldsBndA(path, chrom: string;
                       needed: HashSet[int64]): Table[int64, SlimFields] =
  ## BND variant: END is meaningless (POS+1 sentinel from preproc); just
  ## capture POS and ID, return endP=0 so output uses '.'.
  var vcf: VCF
  if not open(vcf, path):
    raise newException(IOError, "cannot reopen slim BND BCF: " & path)
  var boffData: seq[int32]
  for v in vcf.query(chrom):
    let off = readBoff(v, boffData)
    if off notin needed: continue
    result[off] = (v.POS, 0'i64, $v.ID)
  vcf.close()

proc resolveIntervalPairs(job: MatchJob; pairs: seq[MatchPair]): seq[MatchResult] =
  if pairs.len == 0: return
  var needA, needB: HashSet[int64]
  for p in pairs:
    needA.incl(p.aOff); needB.incl(p.bOff)
  let fieldsA = collectFieldsA(job.pathA, job.chrom, needA)
  var fieldsB: Table[int64, SlimFields]
  for _, pathB in job.binsB:
    for off, f in collectFieldsA(pathB, job.chrom, needB):
      fieldsB[off] = f
  for p in pairs:
    let fA = fieldsA[p.aOff]
    let fB = fieldsB[p.bOff]
    result.add(MatchResult(
      chromA:     job.chrom,
      posA:       fA.pos, endA: fA.endP, idA: fA.id,
      chromB:     job.chrom,
      posB:       fB.pos, endB: fB.endP, idB: fB.id,
      svtype:     job.svtype,
      similarity: p.sim,
      aOffset:    p.aOff, bOffset: p.bOff,
    ))

proc resolveBndPairs(job: MatchJob; pairs: seq[MatchPair]): seq[MatchResult] =
  if pairs.len == 0: return
  var needA, needB: HashSet[int64]
  for p in pairs:
    needA.incl(p.aOff); needB.incl(p.bOff)
  let fieldsA = collectFieldsBndA(job.pathA, job.chrom, needA)
  let fieldsB = collectFieldsBndA(job.binsB[0], job.chrom, needB)
  for p in pairs:
    let fA = fieldsA[p.aOff]
    let fB = fieldsB[p.bOff]
    result.add(MatchResult(
      chromA:     job.chrom,
      posA:       fA.pos, endA: 0, idA: fA.id,
      chromB:     job.chrom,
      posB:       fB.pos, endB: 0, idB: fB.id,
      svtype:     svBND,
      similarity: p.sim,
      aOffset:    p.aOff, bOffset: p.bOff,
    ))

# ---------------------------------------------------------------------------
# Per-job entrypoints (used by the thread pool and direct callers)
# ---------------------------------------------------------------------------

proc runMatchJob*(job: MatchJob, cfg: MatchConfig): seq[MatchResult] =
  ## Stream pairs from matchcore, then resolve fields from the slim BCFs.
  ## Self-mode dedup happens inside matchcore.
  let pairs = streamJobPairs(job, cfg)
  resolveIntervalPairs(job, pairs)

proc runBndMatchJob*(job: MatchJob, cfg: MatchConfig): seq[MatchResult] =
  let pairs = streamBndJobPairs(job, cfg)
  resolveBndPairs(job, pairs)

proc dispatchMatchJob*(job: MatchJob, cfg: MatchConfig): seq[MatchResult] {.inline.} =
  if job.svtype == svBND: runBndMatchJob(job, cfg)
  else:                   runMatchJob(job, cfg)

proc dispatchPairJob*(job: MatchJob, cfg: MatchConfig): seq[MatchPair] {.inline.} =
  ## Pair-only dispatch: callers that don't need resolved fields (collapse).
  if job.svtype == svBND: streamBndJobPairs(job, cfg)
  else:                   streamJobPairs(job, cfg)

# ---------------------------------------------------------------------------
# Thread pool — full MatchResult (match-mode TSV)
# ---------------------------------------------------------------------------

var gJobs:        seq[MatchJob]
var gCfg:         MatchConfig
var gNextJob:     Atomic[int]
var gJobResults:  seq[seq[MatchResult]]

proc threadWorker(dummy: int) {.thread.} =
  {.cast(gcsafe).}:
    while true:
      let idx = gNextJob.fetchAdd(1, moRelaxed)
      if idx >= gJobs.len:
        break
      gJobResults[idx] = dispatchMatchJob(gJobs[idx], gCfg)
      logV("job " & gJobs[idx].chrom & "/" & $gJobs[idx].svtype &
           "/bin" & $gJobs[idx].binA &
           ": " & $gJobResults[idx].len & " matches")

proc runMatchJobsWithPool*(jobs: seq[MatchJob]; cfg: MatchConfig): seq[seq[MatchResult]] =
  ## Run matching jobs via the global thread pool (or inline for nThreads=1).
  result = newSeq[seq[MatchResult]](jobs.len)
  if jobs.len == 0: return
  if cfg.nThreads == 1:
    for i, job in jobs.pairs:
      result[i] = dispatchMatchJob(job, cfg)
      logV("job " & job.chrom & "/" & $job.svtype & "/bin" & $job.binA &
           ": " & $result[i].len & " matches")
  else:
    logV("starting " & $cfg.nThreads & " worker thread(s) for " & $jobs.len & " job(s)")
    gJobs       = jobs
    gCfg        = cfg
    gJobResults = newSeq[seq[MatchResult]](jobs.len)
    gNextJob.store(0, moRelaxed)
    var threads = newSeq[Thread[int]](cfg.nThreads)
    for i in 0 ..< cfg.nThreads:
      createThread(threads[i], threadWorker, i)
    for i in 0 ..< cfg.nThreads:
      joinThread(threads[i])
    result = gJobResults
    logV("workers complete")

# ---------------------------------------------------------------------------
# Thread pool — pair-only (collapse)
# ---------------------------------------------------------------------------

var gPairJobs:     seq[MatchJob]
var gPairCfg:      MatchConfig
var gPairNext:     Atomic[int]
var gPairResults:  seq[seq[MatchPair]]

proc pairThreadWorker(dummy: int) {.thread.} =
  {.cast(gcsafe).}:
    while true:
      let idx = gPairNext.fetchAdd(1, moRelaxed)
      if idx >= gPairJobs.len: break
      gPairResults[idx] = dispatchPairJob(gPairJobs[idx], gPairCfg)
      logV("job " & gPairJobs[idx].chrom & "/" & $gPairJobs[idx].svtype &
           "/bin" & $gPairJobs[idx].binA &
           ": " & $gPairResults[idx].len & " pairs")

proc runMatchPairJobsWithPool*(jobs: seq[MatchJob]; cfg: MatchConfig): seq[seq[MatchPair]] =
  ## Like runMatchJobsWithPool but returns minimal MatchPair triples without
  ## the slim-BCF resolution pass. Used by collapse, which only consumes
  ## (aOff, bOff, sim) via buildSimilarityMap.
  result = newSeq[seq[MatchPair]](jobs.len)
  if jobs.len == 0: return
  if cfg.nThreads == 1:
    for i, job in jobs.pairs:
      result[i] = dispatchPairJob(job, cfg)
      logV("job " & job.chrom & "/" & $job.svtype & "/bin" & $job.binA &
           ": " & $result[i].len & " pairs")
  else:
    logV("starting " & $cfg.nThreads & " pair worker thread(s) for " & $jobs.len & " job(s)")
    gPairJobs    = jobs
    gPairCfg     = cfg
    gPairResults = newSeq[seq[MatchPair]](jobs.len)
    gPairNext.store(0, moRelaxed)
    var threads = newSeq[Thread[int]](cfg.nThreads)
    for i in 0 ..< cfg.nThreads:
      createThread(threads[i], pairThreadWorker, i)
    for i in 0 ..< cfg.nThreads:
      joinThread(threads[i])
    result = gPairResults
    logV("pair workers complete")

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

  let jobResults = runMatchJobsWithPool(jobs, cfg)

  # Write the header line plus all matches in deterministic (job-sorted) order.
  let outFile =
    if cfg.outputPath == "": stdout
    else: open(cfg.outputPath, fmWrite)

  outFile.writeLine("##matcha_metric=" & $cfg.metric)
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
