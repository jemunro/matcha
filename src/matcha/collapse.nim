## collapse.nim — `matcha collapse` subcommand.
##
## Single-pass pipeline (fused preproc + merge via synced_bcf_reader):
##   1. resolveHeaders — analyse N input headers, produce MergedHeader.
##   2. buildFinalHdr  — build shared output header (collapse-specific provenance).
##   3. integratedMerge (mergecore) — stream all N caller VCFs in lockstep via
##      one bcf_srs_t, normalize + filter + write per-(svtype, bin) merged BCFs.
##      Each record gets SRC_INDEX (global sequential) and CALLER_IDX (caller idx).
##   4. runMatchPairJobsWithPool (self-mode, emitSingletons=true) — Pass 1.
##   5. Build allOffsets from MatchPairs (singletons included).
##   6. Scan merged BCFs for cluster members → passQualMap (replaces exploreMerged).
##   7. selfMatchAndCluster → selectRepresentative — cluster and pick representatives.
##   8. writeOutput — stream merged BCFs, buffer representative records,
##      apply output-time INFO filter, sort by coordinate, write final VCF/BCF.

import std/[algorithm, os, sequtils, sets, strutils, tables]
import hts
import hts/private/hts_concat
import utils, preproc, matchcore, mergecore, log, synced_bcf_reader

# ---------------------------------------------------------------------------
# CollapseConfig
# ---------------------------------------------------------------------------

type
  CollapseConfig* = object
    metric*:       Metric
    threshold*:    float64
    bndSlop*:      int
    insSlop*:      int
    insMinSim*:    float64
    linkage*:      LinkageMethod
    priority*:     seq[PriorityCriterion]
    formatFields*: seq[string]   ## FORMAT fields to carry; default ["GT"]
    infoFields*:   seq[string]   ## --info filter; empty = keep only auto-extracted fields
    outputPath*:   string
    nThreads*:     int
    tmpDir*:       string
    callers*:      seq[CallerInput]
    keptChrs*:     seq[string]    ## --chrs: active set; empty = all input contigs.
    chrSet*:       seq[string]    ## --chr-set: universe; empty = all input contigs.
    chunkSize*:    int64          ## --chunk-size: A-side POS range per job.
    writeIndex*:   bool           ## --write-index: emit CSI alongside output.

# ---------------------------------------------------------------------------
# buildFinalHdr — collapse-specific shared output header
# ---------------------------------------------------------------------------

proc buildFinalHdr(callers: seq[CallerInput]; mh: MergedHeader;
                   cfg: CollapseConfig; version, cmdLine: string;
                   outSampleName: var string): ptr bcf_hdr_t =
  ## Build ONE shared header used for: (a) all per-(svtype, bin) merged BCFs
  ## written by integratedMerge, (b) the final output VCF/BCF written by
  ## writeOutput. Eliminates `bcf_translate` at writeOutput time.
  result = bcf_hdr_init("w".cstring)

  addContigsUnion(result, callers, toHashSet(cfg.chrSet))
  addFiltersUnion(result, callers)

  # INFO/FORMAT from mh.headerLines, filtered. INFO is kept only when
  # cfg.infoFields is non-empty AND the field matches the user list (or is
  # one of the always-keep matcha-internal/SV fields); standard SV defs are
  # added below via addStandardSvInfoDefs for the empty-list case.
  let fmtKeep = toHashSet(cfg.formatFields)
  let userInfo = cfg.infoFields
  addHeaderLinesFiltered(result, mh,
    keepInfo = proc (id: string): bool =
      userInfo.len > 0 and keepInfoForMerged(id, userInfo),
    keepFmt  = proc (id: string): bool =
      id in fmtKeep)

  addStandardSvInfoDefs(result)

  # Matcha-internal INFO defs (needed for bcf_translate during merge).
  discard bcf_hdr_append(result,
    "##INFO=<ID=SRC_INDEX,Number=1,Type=Integer,Description=\"matcha-internal: sequential record index\">".cstring)
  discard bcf_hdr_append(result,
    "##INFO=<ID=CALLER_IDX,Number=1,Type=Integer,Description=\"matcha-internal: caller index (0-based)\">".cstring)
  # Provenance INFO defs.
  discard bcf_hdr_append(result,
    "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"Caller names in cluster: representative first, then others (CLI order)\">".cstring)
  discard bcf_hdr_append(result,
    "##INFO=<ID=N_CALLERS,Number=1,Type=Integer,Description=\"Distinct input callsets in cluster\">".cstring)
  discard bcf_hdr_append(result,
    "##INFO=<ID=N_MERGED,Number=1,Type=Integer,Description=\"Total records merged into cluster\">".cstring)
  discard bcf_hdr_append(result,
    ("##source=matcha collapse " & version).cstring)
  if cmdLine.len > 0:
    discard bcf_hdr_append(result, ("##matcha_cmdline=" & cmdLine).cstring)

  # Samples: zero if formatFields is empty, else caller 0's first sample
  # (warn if other callers' first sample differs).
  outSampleName = ""
  if cfg.formatFields.len > 0:
    var vcf0: VCF
    if open(vcf0, callers[0].path):
      let nsamp = bcf_hdr_nsamples(vcf0.header.hdr).int
      if nsamp > 0:
        let samplesArr = cast[cstringArray](vcf0.header.hdr.samples)
        outSampleName = $samplesArr[0]
      vcf0.close()
    for ci in 1 ..< callers.len:
      var vci: VCF
      if not open(vci, callers[ci].path): continue
      let nsampI = bcf_hdr_nsamples(vci.header.hdr).int
      if nsampI > 0:
        let s = $cast[cstringArray](vci.header.hdr.samples)[0]
        if s != outSampleName:
          logWarn("collapse: caller '" & callers[ci].name &
                  "' first sample '" & s & "' differs from caller 0 '" &
                  outSampleName & "' — output uses caller 0's name")
      vci.close()
    if outSampleName.len > 0:
      discard bcf_hdr_add_sample(result, outSampleName.cstring)

  # NULL-add-sample triggers sync; do unconditionally to finalise.
  discard bcf_hdr_add_sample(result, cstring(nil))
  discard bcf_hdr_sync(result)

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

type
  BufferedRep = object
    pos: int64
    rec: ptr bcf1_t  ## owned; caller must bcf_destroy

  ClusterProv = object
    callersStr: string  # comma-joined names: representative first, then others
    nCallers:   int
    nMerged:    int

  CollapseWriterState = object
    outVcf:       VCF
    totalWritten: int

proc keepInfoOut(name: string; infoFilter: seq[string]): bool =
  ## Output-time INFO filter: SRC_INDEX and CALLER_IDX are internal — always drop.
  if name in ["SRC_INDEX", "CALLER_IDX"]: return false
  if name in ["CALLERS", "N_CALLERS", "N_MERGED"]: return true
  if infoFilter.len == 0:
    return name in ["SVTYPE", "SVLEN", "END", "CHR2", "POS2"]
  for tok in infoFilter:
    if name == tok or name.startsWith(tok & "_"): return true
  false

proc openCollapseWriter(cfg: CollapseConfig;
                        finalHdr: ptr bcf_hdr_t;
                        pool: var htsThreadPool): CollapseWriterState =
  ## Open the output VCF/BCF, attach the shared BGZF compression pool,
  ## install the shared finalHdr, and write the header.
  let outPath = if isStdoutPath(cfg.outputPath): "/dev/stdout" else: cfg.outputPath
  let mode =
    if cfg.outputPath.endsWith(".bcf"):      "wb"
    elif cfg.outputPath.endsWith(".vcf.gz"): "wz"
    else:                                    "w"

  if not open(result.outVcf, outPath, mode = mode):
    raise newException(IOError, "cannot open output: " & outPath)
  attachWriterPool(result.outVcf, pool)
  block:
    var dummy: VCF
    discard open(dummy, cfg.callers[0].path)
    result.outVcf.copy_header(dummy.header)
    dummy.close()
  bcf_hdr_destroy(result.outVcf.header.hdr)
  result.outVcf.header.hdr = finalHdr  # SHARED — clear before close to avoid double-free
  discard result.outVcf.write_header()

proc writeCollapseChrom(state: var CollapseWriterState;
                        cfg: CollapseConfig;
                        finalHdr: ptr bcf_hdr_t;
                        mergedBcfs: var seq[VCF];
                        chromName: string;
                        finalClusters: seq[seq[int32]];
                        passQualMap: Table[int32, tuple[hasPASS: bool; qual: uint16; callerIdx: int32]]) =
  ## Stream this chrom's representatives from the open merged BCFs, set
  ## provenance fields, apply the output-time INFO filter, sort by pos, write.
  ## All per-chrom state goes out of scope at proc exit.
  if finalClusters.len == 0: return

  # Build repProv: representative SRC_INDEX → ClusterProv (this chrom's reps only).
  var repProv: Table[int32, ClusterProv]
  for cl in finalClusters:
    if cl.len == 0: continue
    let repIdx = cl[0]
    let repCallerIdx = passQualMap.getOrDefault(repIdx, (false, 0'u16, 0'i32)).callerIdx
    var callerIdxSeen: seq[int32]
    for idx in cl:
      let ci = passQualMap.getOrDefault(idx, (false, 0'u16, 0'i32)).callerIdx
      if ci notin callerIdxSeen: callerIdxSeen.add(ci)
    var callers: seq[string]
    callers.add(cfg.callers[repCallerIdx].name)
    var others: seq[int32]
    for ci in callerIdxSeen:
      if ci != repCallerIdx: others.add(ci)
    others.sort()
    for ci in others: callers.add(cfg.callers[ci].name)
    repProv[repIdx] = ClusterProv(
      callersStr: callers.join(","),
      nCallers:   callerIdxSeen.len,
      nMerged:    cl.len,
    )

  # Query each merged BCF for this chrom, harvest representatives.
  var buf: seq[BufferedRep]
  var idxData: seq[int32]
  for vcf in mergedBcfs.mitems:
    for v in vcf.query(chromName):
      let si = readSrcIndex(v, idxData)
      if si notin repProv: continue
      let prov = repProv[si]

      var callersStr = prov.callersStr
      discard v.info.set("CALLERS", callersStr)
      var nCallers = prov.nCallers.int32
      discard v.info.set("N_CALLERS", nCallers)
      var nMrg = prov.nMerged.int32
      discard v.info.set("N_MERGED", nMrg)

      var toDel: seq[string]
      for fld in v.info.fields:
        if not keepInfoOut(fld.name, cfg.infoFields): toDel.add(fld.name)
      for name in toDel:
        discard v.info.delete(name)

      # Field-ID remap from merged BCF's parsed header to finalHdr (usually
      # a no-op since merged BCFs were written from finalHdr, but the parsed
      # header is a fresh object so be defensive).
      discard bcf_translate(finalHdr, vcf.header.hdr, v.c)

      buf.add(BufferedRep(pos: v.POS, rec: bcf_dup(v.c)))

  # Sort by pos (single chrom).
  buf.sort(proc(a, b: BufferedRep): int = cmp(a.pos, b.pos))

  for br in buf:
    discard bcf_write(vcfHtsFile(state.outVcf), finalHdr, br.rec)
    bcf_destroy(br.rec)
  state.totalWritten += buf.len

proc closeCollapseWriter(state: var CollapseWriterState; cfg: CollapseConfig) =
  ## Close the writer, log, and (if --write-index) build a CSI index.
  state.outVcf.header.hdr = nil   # shared finalHdr; caller owns it
  state.outVcf.close()
  logInfo("collapse: wrote " & $state.totalWritten & " record(s) to " &
          (if isStdoutPath(cfg.outputPath): "stdout" else: cfg.outputPath))
  if cfg.writeIndex and not isStdoutPath(cfg.outputPath) and
     (cfg.outputPath.endsWith(".bcf") or cfg.outputPath.endsWith(".vcf.gz")):
    bcfBuildIndex(cfg.outputPath, cfg.outputPath & ".csi",
                  csi = true, threads = max(cfg.nThreads, 1))
    logInfo("indexed " & cfg.outputPath)

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

proc runCollapse*(cfg: CollapseConfig; cmdLine: string = "") =
  logInfo("matcha collapse: " & $cfg.callers.len & " caller(s)" &
          " linkage=" & $cfg.linkage & " threads=" & $cfg.nThreads)

  for caller in cfg.callers:
    if not fileExists(caller.path):
      logError("input file not found: " & caller.path)
      quit(1)
    if not fileExists(caller.path & ".csi") and not fileExists(caller.path & ".tbi"):
      logError("no index found for: " & caller.path &
               " (run: bcftools index " & caller.path & ")")
      quit(1)
  if cfg.tmpDir == "":
    logError("tmpDir must be set")
    quit(1)

  # Sample-count validation: collapse assumes single-sample inputs (one
  # biological sample per caller). Reject >1 sample columns or inconsistent
  # sample counts across callers (those reflect a cohort, not the single-
  # sample case collapse is designed for).
  block:
    var counts: seq[int]
    for caller in cfg.callers:
      var vcf: VCF
      if not open(vcf, caller.path):
        logError("cannot open: " & caller.path)
        quit(1)
      counts.add(bcf_hdr_nsamples(vcf.header.hdr).int)
      vcf.close()
    for i, n in counts.pairs:
      if n > 1:
        logError("caller '" & cfg.callers[i].name & "' (" &
                 cfg.callers[i].path & ") has " & $n &
                 " sample columns; matcha collapse supports at most " &
                 "1 sample per input (split multi-sample VCFs first)")
        quit(1)
    let first = counts[0]
    for i, n in counts.pairs:
      if n != first:
        logError("inconsistent sample counts across inputs: " &
                 "caller '" & cfg.callers[0].name & "' has " & $first &
                 " sample(s) but caller '" & cfg.callers[i].name &
                 "' has " & $n &
                 "; all collapse inputs must have the same sample count")
        quit(1)

  # Phase 1: resolve output header from all N input headers.
  logVerbose("resolving headers")
  let mh = resolveHeaders(cfg.callers,
    infoFilter = cfg.infoFields,
    fmtFilter  = cfg.formatFields)
  for w in mh.warnings: logWarn("collapse header: " & w)

  # First caller's chrom order anchors output ordering.
  var orderVcf: VCF
  if not open(orderVcf, cfg.callers[0].path):
    raise newException(IOError, "cannot open: " & cfg.callers[0].path)
  let chromOrder = captureChromOrder(orderVcf.header)
  orderVcf.close()

  # Phase 2: build shared output header, then integrated preproc+merge.
  var outSampleName: string
  let finalHdr = buildFinalHdr(cfg.callers, mh, cfg, MatchaVersion, cmdLine,
                               outSampleName)
  logInfo("integrated preproc+merge over " & $cfg.callers.len & " caller(s)")
  let msc = MergeStreamConfig(formatFields:   cfg.formatFields,
                               nThreads:       cfg.nThreads,
                               tmpDir:         cfg.tmpDir,
                               preserveBndAlt: true,
                               preserveInsAlt: true,
                               keptChrs:       toHashSet(cfg.keptChrs),
                               chrSet:         toHashSet(cfg.chrSet))
  let im = integratedMerge(cfg.callers, mh, finalHdr, msc, chromOrder)
  if cfg.keptChrs.len > 0:
    var seen: HashSet[string]
    for caller in cfg.callers:
      var v: VCF
      if not open(v, caller.path): continue
      for c in captureChromOrder(v.header): seen.incl(c)
      v.close()
    warnMissingChrs(cfg.keptChrs, seen)

  # Build a PreprocOutput describing the merged slim BCFs for buildWorkQueue.
  var mergedChromLens: Table[string, int64]
  if im.paths.len > 0:
    var v: VCF
    if open(v, im.paths.values.toSeq[0]):
      for ctg in v.contigs: mergedChromLens[ctg.name] = ctg.length
      v.close()
  let mergedPreproc = PreprocOutput(
    paths:      im.paths,
    populated:  im.populated,
    chromOrder: chromOrder,
    chromLens:  mergedChromLens,
  )

  # Phase 3+4: per-chrom self-match → cluster → emit. Each chrom's clusters
  # are written to the output and dropped before the next chrom is processed,
  # so peak memory tracks the largest single chrom.
  logInfo("self-matching merged slim BCFs")
  let matchCfg = MatchConfig(
    metric:         cfg.metric,
    threshold:      cfg.threshold,
    bndSlop:        cfg.bndSlop,
    insSlop:        cfg.insSlop,
    insMinSim:      cfg.insMinSim,
    nThreads:       cfg.nThreads,
    tmpDir:         cfg.tmpDir,
    selfMode:       true,
    emitSingletons: true,
    chunkSize:      cfg.chunkSize,
  )

  let (allJobs, _) = buildWorkQueue(mergedPreproc, mergedPreproc, matchCfg)
  var jobsByChrom: Table[int16, seq[MatchJob]]
  for j in allJobs:
    jobsByChrom.mgetOrPut(j.chromIdx, @[]).add(j)

  var writerPool = newWriterPool(cfg.nThreads)
  var writerState = openCollapseWriter(cfg, finalHdr, writerPool)

  # Open merged slim BCFs once; each is CSI-indexed by integratedMerge.
  var mergedBcfs: seq[VCF]
  for path in im.paths.values:
    var v: VCF
    if not open(v, path):
      raise newException(IOError, "cannot open merged BCF: " & path)
    mergedBcfs.add(v)

  for chromIdx in 0'i16 ..< chromOrder.len.int16:
    let cjobs = jobsByChrom.getOrDefault(chromIdx, @[])
    if cjobs.len == 0: continue
    var locByIdx: Table[int32, tuple[chromIdx: int16; pos: int32; fileIdx: int16]]
    var passQualMap: Table[int32, tuple[hasPASS: bool; qual: uint16; callerIdx: int32]]
    var finalClusters: seq[seq[int32]]
    selfMatchAndClusterChrom(
      mergedPreproc, matchCfg, cjobs,
      cfg.linkage, cfg.threshold, cfg.priority,
      "collapse self-match", warnCallerStats = true,
      locByIdx, passQualMap, finalClusters)
    writeCollapseChrom(writerState, cfg, finalHdr, mergedBcfs,
                       chromOrder[chromIdx], finalClusters, passQualMap)
    # locByIdx, passQualMap, finalClusters go out of scope → freed.

  for v in mergedBcfs.mitems: v.close()
  closeCollapseWriter(writerState, cfg)
  destroyWriterPool(writerPool)

  # Clean up: per-invocation temp dir (merged slim BCFs + CSI indexes), and finalHdr.
  removeDir(cfg.tmpDir)
  bcf_hdr_destroy(finalHdr)
