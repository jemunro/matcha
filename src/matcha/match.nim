## match.nim — per-job matching logic, thread pool, and output assembly.

import std/[algorithm, atomics, os, strutils, tables]
import hts
import utils, intervals, preproc, log

# ---------------------------------------------------------------------------
# Per-job matching
# ---------------------------------------------------------------------------

proc processJob*(job: MatchJob, cfg: MatchConfig): int {.discardable.} =
  ## Stream A records for this (chrom, svtype), query the B index for
  ## candidates, compute metrics, and write matching rows to job.outTsv.
  ## Returns the number of matches written.
  var vcfA, vcfB: VCF
  if not open(vcfA, job.pathA):
    raise newException(IOError, "cannot open A BCF: " & job.pathA)
  if not open(vcfB, job.pathB):
    raise newException(IOError, "cannot open B BCF: " & job.pathB)

  var outFile = open(job.outTsv, fmWrite)

  let effectiveThreshold =
    if cfg.minOverlapSet and cfg.minJaccardSet: max(cfg.minOverlap, cfg.minJaccard)
    elif cfg.minOverlapSet: cfg.minOverlap
    else: cfg.minJaccard

  var svlenData: seq[int32]
  var endData:   seq[int32]

  for va in vcfA.query(job.chrom):
    let hasEnd  = va.info().get("END",   endData)   == Status.OK and endData.len   > 0
    let hasSvln = va.info().get("SVLEN", svlenData) == Status.OK and svlenData.len > 0

    var svlenA: int64
    if hasSvln:
      svlenA = abs(int64(svlenData[0]))
    elif hasEnd:
      svlenA = int64(endData[0]) - va.POS
    else:
      continue

    let endA = if hasEnd: int64(endData[0]) else: va.POS + svlenA
    let win  = queryWindow(svlenA, effectiveThreshold)
    let regionStr = $va.CHROM & ":" & $max(1'i64, va.POS - win) & "-" & $(endA + win)

    for vb in vcfB.query(regionStr):
      let hasEndB  = vb.info().get("END",   endData)   == Status.OK and endData.len   > 0
      let hasSvlnB = vb.info().get("SVLEN", svlenData) == Status.OK and svlenData.len > 0

      var svlenB: int64
      if hasSvlnB:
        svlenB = abs(int64(svlenData[0]))
      elif hasEndB:
        svlenB = int64(endData[0]) - vb.POS
      else:
        continue

      let endB = if hasEndB: int64(endData[0]) else: vb.POS + svlenB

      let ovl = reciprocalOverlap(va.POS, endA, vb.POS, endB)
      let jac = jaccard(va.POS, endA, vb.POS, endB)

      let passOverlap = (not cfg.minOverlapSet) or (ovl >= cfg.minOverlap)
      let passJaccard = (not cfg.minJaccardSet) or (jac >= cfg.minJaccard)

      if passOverlap and passJaccard:
        outFile.writeLine(formatMatchResult(MatchResult(
          chrom:    $va.CHROM,
          posA:     va.POS,
          endA:     endA,
          idA:      $va.ID,
          posB:     vb.POS,
          endB:     endB,
          idB:      $vb.ID,
          svtype:   job.svtype,
          overlap:  ovl,
          jaccard:  jac,
        )))
        inc result

  outFile.close()
  vcfA.close()
  vcfB.close()

# ---------------------------------------------------------------------------
# Thread pool (global state + atomic counter)
# ---------------------------------------------------------------------------

var gJobs:    seq[MatchJob]
var gCfg:     MatchConfig
var gNextJob: Atomic[int]

proc threadWorker(dummy: int) {.thread.} =
  {.cast(gcsafe).}:
    while true:
      let idx = gNextJob.fetchAdd(1, moRelaxed)
      if idx >= gJobs.len:
        break
      let n = processJob(gJobs[idx], gCfg)
      logV("job " & gJobs[idx].chrom & "/" & $gJobs[idx].svtype &
           ": " & $n & " matches")

# ---------------------------------------------------------------------------
# Parallel preprocessing (used when nThreads >= 2)
# ---------------------------------------------------------------------------
# A and B preprocessing are independent: they read different inputs and write
# to different temp paths (distinct "A"/"B" prefixes). Two threads, each
# running preprocessVcf, are joined before the work queue is built.

var gPpInputs:  array[2, tuple[path, tmpDir, prefix: string]]
var gPpOutputs: array[2, PreprocOutput]

proc preprocWorker(idx: int) {.thread.} =
  {.cast(gcsafe).}:
    let s = gPpInputs[idx]
    gPpOutputs[idx] = preprocessVcf(s.path, s.tmpDir, s.prefix)

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

proc runMatch*(cfg: MatchConfig) =
  logV("matcha match: A=" & cfg.callsetA & " B=" & cfg.callsetB &
       " threads=" & $cfg.nThreads & " tmp=" & cfg.tmpDir)

  var filesA, filesB: PreprocOutput
  if cfg.nThreads >= 2:
    logV("preprocessing A and B in parallel")
    gPpInputs[0] = (cfg.callsetA, cfg.tmpDir, "A")
    gPpInputs[1] = (cfg.callsetB, cfg.tmpDir, "B")
    var thA, thB: Thread[int]
    createThread(thA, preprocWorker, 0)
    createThread(thB, preprocWorker, 1)
    joinThread(thA)
    joinThread(thB)
    filesA = gPpOutputs[0]
    filesB = gPpOutputs[1]
  else:
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A")
    filesB = preprocessVcf(cfg.callsetB, cfg.tmpDir, "B")

  var jobs = buildWorkQueue(filesA, filesB, cfg.tmpDir)
  logV("work queue: " & $jobs.len & " (chrom,svtype) job(s)")
  if jobs.len == 0:
    return

  # Sort for deterministic output order (chrom then svtype string)
  jobs.sort do (a, b: MatchJob) -> int:
    if a.chrom != b.chrom: cmp(a.chrom, b.chrom)
    else: cmp($a.svtype, $b.svtype)

  if cfg.nThreads == 1:
    for job in jobs:
      let n = processJob(job, cfg)
      logV("job " & job.chrom & "/" & $job.svtype & ": " & $n & " matches")
  else:
    logV("starting " & $cfg.nThreads & " worker thread(s)")
    gJobs = jobs
    gCfg  = cfg
    gNextJob.store(0, moRelaxed)
    var threads = newSeq[Thread[int]](cfg.nThreads)
    for i in 0 ..< cfg.nThreads:
      createThread(threads[i], threadWorker, i)
    for i in 0 ..< cfg.nThreads:
      joinThread(threads[i])
    logV("workers complete")

  # Concatenate per-job TSVs in sorted order, prefixed with a header line.
  let outFile =
    if cfg.outputPath == "": stdout
    else: open(cfg.outputPath, fmWrite)

  outFile.writeLine(OutputHeader)

  var totalMatches = 0
  for job in jobs:
    if fileExists(job.outTsv):
      let content = readFile(job.outTsv)
      if content.len > 0:
        outFile.write(content)
        totalMatches += content.count('\n')
      removeFile(job.outTsv)

  if cfg.outputPath != "":
    outFile.close()
  logV("wrote " & $totalMatches & " match(es)" &
       (if cfg.outputPath != "": " to " & cfg.outputPath else: " to stdout"))

  # Clean up temp BCFs
  for path in filesA.paths.values:
    if fileExists(path):     removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")
  for path in filesB.paths.values:
    if fileExists(path):     removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")
