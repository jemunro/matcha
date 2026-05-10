## match.nim — per-job matching with size bins + tiled B buffers, thread pool,
## and output assembly.

import std/[algorithm, atomics, os, sequtils, tables]
import hts
import utils, intervals, preproc, log, bins

# ---------------------------------------------------------------------------
# Per-job matching
# ---------------------------------------------------------------------------

proc readBoff(v: Variant, scratch: var seq[int32]): int64 =
  ## Decode INFO/MATCHA_BOFF (Number=2 Integer: high32, low32) into an int64.
  ## Returns 0 if the field is absent. The mask on the low half avoids sign-
  ## extension when the int32 has its high bit set.
  if v.info().get("MATCHA_BOFF", scratch) != Status.OK or scratch.len < 2:
    return 0
  result = (int64(scratch[0]) shl 32) or (int64(scratch[1]) and 0xFFFFFFFF'i64)

proc extractEnd(v: Variant, endData, svlenData: var seq[int32];
                outEnd: var int64): bool =
  ## Resolve END for a slim record. After preproc, END is always written
  ## authoritatively, so this is a fast path; SVLEN fallback is defensive.
  if v.info().get("END", endData) == Status.OK and endData.len > 0:
    outEnd = int64(endData[0]); return true
  if v.info().get("SVLEN", svlenData) == Status.OK and svlenData.len > 0:
    outEnd = v.POS + abs(int64(svlenData[0])); return true
  false

proc runMatchJob*(job: MatchJob, cfg: MatchConfig): seq[MatchResult] =
  ## Stream A records from the per-(svtype, binA) BCF restricted to job.chrom.
  ## For each A record, query each adjacent populated B bin via the per-binB
  ## TiledBuffer. Apply the metric filter; emit MatchResults carrying both
  ## sides' MATCHA_BOFF source-file offsets for downstream modes.
  var vcfA: VCF
  if not open(vcfA, job.pathA):
    raise newException(IOError, "cannot open A BCF: " & job.pathA)

  # Lazily-opened per-binB readers + tiled buffers. Sorted bin order so the
  # output is deterministic regardless of hash-table iteration order.
  var vcfsB: Table[int, VCF]
  var buffers: Table[int, TiledBuffer]
  var sortedBinsB = toSeq(job.binsB.keys)
  sortedBinsB.sort()
  for binB in sortedBinsB:
    buffers[binB] = initTiledBuffer(binRange(binB).hi, job.chrom)

  # Scratch buffers for INFO field decodes.
  var endData, svlenData, boffData: seq[int32]

  # Region-query one tile from the per-binB slim BCF. Lazily opens the reader
  # on first call for each binB. CSI queries return records whose [POS, END)
  # overlaps the region, so a record straddling a tile boundary appears in both
  # adjacent fetches; filter by POS to assign each record to exactly one tile.
  proc fetchTile(binB, tileIdx: int): seq[BufferedRec] =
    if binB notin vcfsB:
      var v: VCF
      if not open(v, job.binsB[binB]):
        raise newException(IOError, "cannot open B BCF: " & job.binsB[binB])
      vcfsB[binB] = v
    let W = binRange(binB).hi
    let regStart = max(1'i64, tileIdx.int64 * W)
    let regEnd   = (tileIdx.int64 + 1) * W - 1
    let region   = job.chrom & ":" & $regStart & "-" & $regEnd
    for vb in vcfsB[binB].query(region):
      if int(vb.POS div W) != tileIdx: continue
      var endB: int64
      if not extractEnd(vb, endData, svlenData, endB): continue
      result.add(BufferedRec(
        pos: vb.POS, endPos: endB,
        id: $vb.ID,
        bOffset: readBoff(vb, boffData),
      ))

  for va in vcfA.query(job.chrom):
    var endA: int64
    if not extractEnd(va, endData, svlenData, endA): continue
    let posA = va.POS
    let aOff = readBoff(va, boffData)

    for binB in sortedBinsB:
      # Position window [posA - U, posA + svlenA = endA). Asymmetric: a B
      # record up to U bp left of posA can extend rightward into A, while
      # B records at posA + svlenA or later cannot overlap.
      let b = binB  # owned copy — lent iterator vars can't be captured in closures
      let cands = buffers[b].getCandidates(posA, endA,
        proc(ti: int): seq[BufferedRec] = fetchTile(b, ti))
      for cand in cands:
        let ovl = reciprocalOverlap(posA, endA, cand.pos, cand.endPos)
        let jac = jaccard(posA, endA, cand.pos, cand.endPos)
        let passOverlap = (not cfg.minOverlapSet) or (ovl >= cfg.minOverlap)
        let passJaccard = (not cfg.minJaccardSet) or (jac >= cfg.minJaccard)
        if passOverlap and passJaccard:
          result.add(MatchResult(
            chrom:   $va.CHROM,
            posA:    posA,
            endA:    endA,
            idA:     $va.ID,
            posB:    cand.pos,
            endB:    cand.endPos,
            idB:     cand.id,
            svtype:  job.svtype,
            overlap: ovl,
            jaccard: jac,
            aOffset: aOff,
            bOffset: cand.bOffset,
          ))

    # Evict tiles no future A record can need (A is position-sorted within
    # this chrom-restricted stream).
    for binB, buf in buffers.mpairs:
      buf.evict(posA)

  vcfA.close()
  for v in vcfsB.mvalues:
    v.close()

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
  for path in filesA.paths.values:
    if fileExists(path):     removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")
  for path in filesB.paths.values:
    if fileExists(path):     removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")
