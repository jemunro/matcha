## setops.nim — `matcha intersect` and `matcha setdiff` subcommands.
##
## Record-level set operations on two SV callsets, sharing the binning +
## tiled-buffer matching engine with `match`/`anno`. Both re-emit A records
## verbatim (genotypes/FORMAT preserved); they differ only in the predicate:
##
##   intersect : keep an A record iff it matches >=1 record in B.
##   setdiff   : keep an A record iff it matches no record in B.
##
## Three phases:
##   1. Preproc A and B (default keep-set; no extra INFO fields).
##   2. Match A vs B, collecting the SRC_INDEX of every A record with a match.
##   3. Stream the *original* A file, re-emitting kept records unchanged.

import std/[os, sets, strutils]
import hts
import hts/private/hts_concat
import utils, preproc, log, match

type SetOpConfig* = object
  metric*:      Metric       ## Active interval metric (mOverlap | mJaccard).
  threshold*:   float64      ## Minimum score for the active metric.
  bndSlop*:     int          ## --bnd-slop (default 50)
  insSlop*:     int          ## --ins-slop (default 50)
  insMinSim*:   float64      ## --min-ins-sim (default 0.75)
  nThreads*:    int
  tmpDir*:      string
  outputPath*:  string
  callsetA*:    string       ## input being filtered
  callsetB*:    string       ## reference set
  keptChrs*:    seq[string]  ## --chrs: active set (output); empty = all input contigs.
  chrSet*:      seq[string]  ## --chr-set: universe; empty = all input contigs.
  chunkSize*:   int64        ## --chunk-size: A-side POS range per job.
  writeIndex*:  bool         ## --write-index: emit CSI alongside output.
  keepMatched*: bool         ## true = intersect, false = setdiff.

proc opName(cfg: SetOpConfig): string =
  if cfg.keepMatched: "intersect" else: "setdiff"

proc validateOutputPath(p: string) =
  if isStdoutPath(p): return
  if not (p.endsWith(".vcf") or p.endsWith(".vcf.gz") or p.endsWith(".bcf")):
    logError("-o must end in .vcf, .vcf.gz, or .bcf: " & p)
    quit(1)

proc setOpMatchCfg(cfg: SetOpConfig): MatchConfig {.inline.} =
  MatchConfig(metric: cfg.metric, threshold: cfg.threshold,
              bndSlop: cfg.bndSlop, insSlop: cfg.insSlop,
              insMinSim: cfg.insMinSim, nThreads: cfg.nThreads,
              tmpDir: cfg.tmpDir, selfMode: false, emitSingletons: false,
              chunkSize: cfg.chunkSize)

proc runSetOp*(cfg: var SetOpConfig) =
  logInfo("matcha " & opName(cfg) & ": A=" & cfg.callsetA & " B=" & cfg.callsetB &
          " threads=" & $cfg.nThreads & " tmp=" & cfg.tmpDir)
  validateOutputPath(cfg.outputPath)

  let keptChrsSet = toHashSet(cfg.keptChrs)
  let chrSetSet   = toHashSet(cfg.chrSet)

  # --- Phase 1: preproc -----------------------------------------------------
  var filesA, filesB: PreprocOutput
  if cfg.nThreads >= 2:
    logInfo("preprocessing A and B in parallel")
    (filesA, filesB) = runParallelPreproc(
      PreprocInput(path: cfg.callsetA, tmpDir: cfg.tmpDir, prefix: "A",
                   ioThreads: 2, keptChrs: keptChrsSet, chrSet: chrSetSet),
      PreprocInput(path: cfg.callsetB, tmpDir: cfg.tmpDir, prefix: "B",
                   ioThreads: 2, keptChrs: keptChrsSet, chrSet: chrSetSet))
  else:
    filesA = preprocessVcf(cfg.callsetA, cfg.tmpDir, "A",
                           keptChrs = keptChrsSet, chrSet = chrSetSet)
    filesB = preprocessVcf(cfg.callsetB, cfg.tmpDir, "B",
                           keptChrs = keptChrsSet, chrSet = chrSetSet)
  var chrsSeen = filesA.chrsSeen
  for c in filesB.chrsSeen: chrsSeen.incl(c)
  warnMissingChrs(cfg.keptChrs, chrsSeen)

  let matchCfg = setOpMatchCfg(cfg)
  let (jobs, _) = buildWorkQueue(filesA, filesB, matchCfg)
  logInfo(opName(cfg) & " work queue: " & $jobs.len & " (chrom, svtype, binA) job(s)")

  # --- Phase 2: match, collect matched A indices ----------------------------
  var matched = initHashSet[int32]()
  if jobs.len > 0:
    let jobResults = runMatchPairJobsWithPool(jobs, matchCfg)
    for jrs in jobResults:
      for p in jrs:
        if p.srcIndexB != NO_MATCH:
          matched.incl(p.srcIndexA)
  logInfo(opName(cfg) & ": " & $matched.len & " A record(s) with a match")

  # --- Phase 3: stream original A, write kept records verbatim --------------
  var vcfA: VCF
  if not open(vcfA, cfg.callsetA, threads = if cfg.nThreads >= 2: 2 else: 0):
    raise newException(IOError, "cannot reopen input: " & cfg.callsetA)
  if chrSetSet.len > 0:
    # Drop contigs outside the universe (chr-set) from the copied output header.
    var n: cint = 0
    let names = bcf_hdr_seqnames(vcfA.header.hdr, n.addr)
    if names != nil:
      var toRemove: seq[string]
      for i in 0 ..< n.int:
        let c = $names[i]
        if c notin chrSetSet: toRemove.add(c)
      free(names)
      for c in toRemove:
        bcf_hdr_remove(vcfA.header.hdr, BCF_HEADER_TYPE.BCF_HL_CTG.cint, c.cstring)
      discard bcf_hdr_sync(vcfA.header.hdr)

  var outWriter: VCF
  # hts-nim's close() skips hts_close() for fname=="-"; route stdout via
  # "/dev/stdout" so the BGZF buffer flushes.
  let outPath = if isStdoutPath(cfg.outputPath): "/dev/stdout" else: cfg.outputPath
  if not open(outWriter, outPath, mode = "w"):
    raise newException(IOError, "cannot open output: " & outPath)
  var writerPool = newWriterPool(cfg.nThreads)
  attachWriterPool(outWriter, writerPool)
  outWriter.copy_header(vcfA.header)
  if not outWriter.write_header():
    raise newException(IOError, "failed to write output header")

  # Stream original A per-chrom (mirroring preproc's CSI iteration), joining on
  # SRC_INDEX — the same sequential counter preproc assigned.
  var srcIdx: int32 = 0
  var nInput = 0
  var nKept = 0
  for chrom in filesA.chromOrder:
    if keptChrsSet.len > 0 and chrom notin keptChrsSet: continue
    for v in vcfA.query(chrom):
      let curIdx = srcIdx
      inc srcIdx
      inc nInput
      let keep = if cfg.keepMatched: curIdx in matched else: curIdx notin matched
      if not keep: continue
      if not outWriter.write_variant(v):
        raise newException(IOError, "failed to write variant at " &
          $v.CHROM & ":" & $v.POS)
      inc nKept

  vcfA.close()
  outWriter.close()
  destroyWriterPool(writerPool)
  logInfo(opName(cfg) & ": wrote " & $nKept & " of " & $nInput & " record(s) to " &
          (if isStdoutPath(cfg.outputPath): "stdout" else: cfg.outputPath))
  if cfg.writeIndex and not isStdoutPath(cfg.outputPath) and
     (cfg.outputPath.endsWith(".bcf") or cfg.outputPath.endsWith(".vcf.gz")):
    bcfBuildIndex(cfg.outputPath, cfg.outputPath & ".csi",
                  csi = true, threads = max(cfg.nThreads, 1))
    logInfo("indexed " & cfg.outputPath)

  # Clean up: per-invocation temp dir (holds A and B preproc artifacts).
  removeDir(cfg.tmpDir)
