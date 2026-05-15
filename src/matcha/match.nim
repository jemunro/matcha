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

proc collectFields(path, chrom: string; needed: HashSet[int64];
                   bnd: bool): Table[int64, SlimFields] =
  ## Scan a slim BCF for `chrom`, capturing (POS, END, ID) keyed by
  ## MATCHA_BOFF for every record whose offset is in `needed`. For BND
  ## records END is the POS+1 sentinel from preproc; return endP=0 so the
  ## output renders as '.'.
  var vcf: VCF
  if not open(vcf, path):
    raise newException(IOError, "cannot reopen slim BCF: " & path)
  var endData, svlenData, boffData: seq[int32]
  for v in vcf.query(chrom):
    let off = readBoff(v, boffData)
    if off notin needed: continue
    if bnd:
      result[off] = (v.POS, 0'i64, $v.ID)
    else:
      var endP: int64
      if not extractEnd(v, endData, svlenData, endP): continue
      result[off] = (v.POS, endP, $v.ID)
  vcf.close()

proc resolveIntervalPairs(job: MatchJob; pairs: seq[MatchPair]): seq[MatchResult] =
  if pairs.len == 0: return
  var needA, needB: HashSet[int64]
  for p in pairs:
    needA.incl(p.aOff); needB.incl(p.bOff)
  let fieldsA = collectFields(job.pathA, job.chrom, needA, bnd = false)
  var fieldsB: Table[int64, SlimFields]
  for _, pathB in job.binsB:
    for off, f in collectFields(pathB, job.chrom, needB, bnd = false):
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
    ))

proc resolveBndPairs(job: MatchJob; pairs: seq[MatchPair]): seq[MatchResult] =
  if pairs.len == 0: return
  var needA, needB: HashSet[int64]
  for p in pairs:
    needA.incl(p.aOff); needB.incl(p.bOff)
  let fieldsA = collectFields(job.pathA, job.chrom, needA, bnd = true)
  let fieldsB = collectFields(job.binsB[0], job.chrom, needB, bnd = true)
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
    ))

# ---------------------------------------------------------------------------
# Per-job entrypoints (used by the thread pool and direct callers)
# ---------------------------------------------------------------------------

proc dispatchPairJob*(job: MatchJob, cfg: MatchConfig): seq[MatchPair] =
  ## Stream (aOff, bOff, sim) triples for one job from matchcore. Used
  ## directly by collapse, and as the first step in dispatchMatchJob.
  ## Self-mode dedup happens inside matchcore.
  if job.svtype == svBND: streamBndJobPairs(job, cfg)
  else:                   streamJobPairs(job, cfg)

proc dispatchMatchJob*(job: MatchJob, cfg: MatchConfig): seq[MatchResult] =
  ## Pairs from dispatchPairJob plus slim-BCF resolution into MatchResults.
  let pairs = dispatchPairJob(job, cfg)
  if job.svtype == svBND: resolveBndPairs(job, pairs)
  else:                   resolveIntervalPairs(job, pairs)

# ---------------------------------------------------------------------------
# Generic per-job thread pool — parameterized on the per-job result type
# ---------------------------------------------------------------------------

type
  JobRunner[R] = proc(job: MatchJob, cfg: MatchConfig): seq[R] {.nimcall.}

  PoolState[R] = object
    jobs:    seq[MatchJob]
    cfg:     MatchConfig
    next:    Atomic[int]
    results: seq[seq[R]]
    runner:  JobRunner[R]
    label:   string   # "matches" / "pairs", for log lines only

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
  ## Dispatch `jobs` across `cfg.nThreads` workers, calling `runner` per job.
  ## Each worker pulls indexes via an atomic counter; results land in the
  ## per-job slot. The state lives on this proc's stack and all threads are
  ## joined before return.
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

proc runMatchJobsWithPool*(jobs: seq[MatchJob]; cfg: MatchConfig): seq[seq[MatchResult]] =
  ## Full MatchResult pipeline: pair streaming + slim-BCF resolution.
  runJobsWithPool[MatchResult](jobs, cfg, dispatchMatchJob, "matches")

proc runMatchPairJobsWithPool*(jobs: seq[MatchJob]; cfg: MatchConfig): seq[seq[MatchPair]] =
  ## Pair-only pipeline: returns (aOff, bOff, sim) triples without the
  ## slim-BCF resolution pass. Used by collapse, which only consumes
  ## triples via buildSimilarityMap.
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
