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
##   7. clusterAll → selectRepresentative — cluster and pick representatives.
##   8. writeOutput — stream merged BCFs, buffer representative records,
##      apply output-time INFO filter, sort by coordinate, write final VCF/BCF.

import std/[algorithm, os, sets, strutils, tables]
import hts
import hts/private/hts_concat
import utils, preproc, match, matchcore, mergecore, log, synced_bcf_reader

# ---------------------------------------------------------------------------
# CollapseConfig
# ---------------------------------------------------------------------------

type
  CollapseConfig* = object
    metric*:       Metric
    threshold*:    float64
    bndSlop*:      int
    linkage*:      LinkageMethod
    priority*:     seq[PriorityCriterion]
    formatFields*: seq[string]   ## FORMAT fields to carry; default ["GT"]
    infoFields*:   seq[string]   ## --info filter; empty = keep only auto-extracted fields
    outputPath*:   string
    nThreads*:     int
    tmpDir*:       string
    callers*:      seq[CallerInput]

# ---------------------------------------------------------------------------
# Header traversal helpers (used by buildFinalHdr)
# ---------------------------------------------------------------------------

proc infoFieldDefs(h: ptr bcf_hdr_t): seq[(string, string, string, string)] =
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_INFO.cint: continue
    let keys = cast[ptr UncheckedArray[cstring]](hr.keys)
    let vals = cast[ptr UncheckedArray[cstring]](hr.vals)
    var id, num, typ, desc = ""
    for j in 0 ..< hr.nkeys.int:
      case $keys[j]
      of "ID":     id  = $vals[j]
      of "Number": num = $vals[j]
      of "Type":   typ = $vals[j]
      of "Description":
        let dv = $vals[j]
        desc = if dv.len >= 2 and dv[0] == '"' and dv[^1] == '"':
                 dv[1 ..< dv.len - 1]
               else: dv
    if id.len > 0: result.add((id, num, typ, desc))

# ---------------------------------------------------------------------------
# Pass 2 — enumerate allOffsets + PASS/QUAL from merged slim BCFs
# ---------------------------------------------------------------------------

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

  # Contigs: union from all callers (first caller wins for order).
  var seenCtg: HashSet[string]
  for caller in callers:
    var vcf: VCF
    if not open(vcf, caller.path): continue
    var n: cint = 0
    let names = bcf_hdr_seqnames(vcf.header.hdr, n.addr)
    if names != nil:
      for i in 0 ..< n.int:
        let c = $names[i]
        if c notin seenCtg:
          seenCtg.incl(c)
          discard bcf_hdr_append(result,
            ("##contig=<ID=" & c & ">").cstring)
      free(names)
    vcf.close()

  # FILTER defs: merge from all callers.
  var seenFlt: HashSet[string]
  for caller in callers:
    var vcf: VCF
    if not open(vcf, caller.path): continue
    for line in collectFilterLines(vcf.header.hdr):
      if line notin seenFlt:
        seenFlt.incl(line)
        discard bcf_hdr_append(result, line.cstring)
    vcf.close()

  # INFO/FORMAT from mh.headerLines, filtered.
  let fmtKeep = toHashSet(cfg.formatFields)
  for line in mh.headerLines:
    let idStart = line.find("ID=")
    if idStart < 0:
      discard bcf_hdr_append(result, line.cstring)
      continue
    let idEnd = line.find(',', idStart + 3)
    let fieldId =
      if idEnd > 0: line[idStart + 3 ..< idEnd]
      else:         line[idStart + 3 ..< line.len - 1]
    if line.startsWith("##FORMAT"):
      if fieldId in fmtKeep:
        discard bcf_hdr_append(result, line.cstring)
    else:
      if cfg.infoFields.len > 0 and keepInfoForMerged(fieldId, cfg.infoFields):
        discard bcf_hdr_append(result, line.cstring)

  # Ensure standard SV INFO defs (in case callers' headers omit any).
  let outHdr = result
  proc hasInfo(name: string): bool =
    for (id, _, _, _) in infoFieldDefs(outHdr):
      if id == name: return true
    false
  if not hasInfo("SVTYPE"):
    discard bcf_hdr_append(result,
      "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">".cstring)
  if not hasInfo("SVLEN"):
    discard bcf_hdr_append(result,
      "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV (absolute value)\">".cstring)
  if not hasInfo("END"):
    discard bcf_hdr_append(result,
      "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV (1-based, inclusive)\">".cstring)
  if not hasInfo("CHR2"):
    discard bcf_hdr_append(result,
      "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome of mate breakend\">".cstring)
  if not hasInfo("POS2"):
    discard bcf_hdr_append(result,
      "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Position of mate breakend\">".cstring)

  # Matcha-internal INFO defs (needed for bcf_translate during merge).
  discard bcf_hdr_append(result,
    "##INFO=<ID=SRC_INDEX,Number=1,Type=Integer,Description=\"matcha-internal: sequential record index\">".cstring)
  discard bcf_hdr_append(result,
    "##INFO=<ID=CALLER_IDX,Number=1,Type=Integer,Description=\"matcha-internal: caller index (0-based)\">".cstring)
  # Provenance INFO defs.
  discard bcf_hdr_append(result,
    "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Representative caller name\">".cstring)
  discard bcf_hdr_append(result,
    "##INFO=<ID=SOURCELIST,Number=.,Type=String,Description=\"All caller names in cluster (CLI order)\">".cstring)
  discard bcf_hdr_append(result,
    "##INFO=<ID=N_SOURCE,Number=1,Type=Integer,Description=\"Distinct input callsets in cluster\">".cstring)
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
    chromOrderIdx: int
    pos:           int64
    rec:           ptr bcf1_t  ## owned; caller must bcf_destroy

  ClusterProv = object
    sourceList: seq[string]
    nSource:    int
    nMerged:    int

proc writeOutput(cfg: CollapseConfig;
                 finalHdr: ptr bcf_hdr_t;
                 chromOrder: seq[string];
                 mergedPaths: Table[SvtypeBin, string];
                 finalClusters: seq[seq[int32]];
                 passQualMap: Table[int32, tuple[hasPASS: bool; qual: float32; callerIdx: int32]]) =
  ## Stream merged BCFs, identify cluster representatives by SRC_INDEX,
  ## set provenance fields, apply output-time INFO filter, sort by coordinate,
  ## write final VCF/BCF.

  # Build repProv: representative SRC_INDEX → ClusterProv.
  var repProv: Table[int32, ClusterProv]
  for cl in finalClusters:
    if cl.len == 0: continue
    let repIdx = cl[0]
    var callerIdxSeen: seq[int32]
    for idx in cl:
      let ci = passQualMap.getOrDefault(idx, (false, 0f32, 0'i32)).callerIdx
      if ci notin callerIdxSeen: callerIdxSeen.add(ci)
    callerIdxSeen.sort()
    var sourceList: seq[string]
    for ci in callerIdxSeen: sourceList.add(cfg.callers[ci].name)
    repProv[repIdx] = ClusterProv(
      sourceList: sourceList,
      nSource:    sourceList.len,
      nMerged:    cl.len,
    )

  # chrom → header-order index for sort.
  var chromIdx: Table[string, int]
  for i, c in chromOrder.pairs: chromIdx[c] = i

  # Output-time INFO filter: SRC_INDEX and CALLER_IDX are internal — always drop.
  let infoFilter = cfg.infoFields
  proc keepInfoOut(name: string): bool =
    if name in ["SRC_INDEX", "CALLER_IDX"]: return false
    if name in ["SOURCE", "SOURCELIST", "N_SOURCE", "N_MERGED"]: return true
    if infoFilter.len == 0:
      return name in ["SVTYPE", "SVLEN", "END", "CHR2", "POS2"]
    for tok in infoFilter:
      if name == tok or name.startsWith(tok & "_"): return true
    false

  # Stream merged BCFs, collect representatives.
  var buf: seq[BufferedRep]
  var idxData: seq[int32]
  for path in mergedPaths.values:
    var vcf: VCF
    if not open(vcf, path):
      raise newException(IOError, "cannot open merged BCF: " & path)

    for v in vcf:
      let si = readSrcIndex(v, idxData)
      if si notin repProv: continue
      let prov = repProv[si]

      var slStr = prov.sourceList.join(",")
      discard v.info.set("SOURCELIST", slStr)
      var nSrc = prov.nSource.int32
      discard v.info.set("N_SOURCE", nSrc)
      var nMrg = prov.nMerged.int32
      discard v.info.set("N_MERGED", nMrg)

      var toDel: seq[string]
      for fld in v.info.fields:
        if not keepInfoOut(fld.name): toDel.add(fld.name)
      for name in toDel:
        discard v.info.delete(name)

      # Field-ID remap from merged BCF's parsed header to finalHdr (usually
      # a no-op since merged BCFs were written from finalHdr, but the parsed
      # header is a fresh object so be defensive).
      discard bcf_translate(finalHdr, vcf.header.hdr, v.c)

      let chromOI = chromIdx.getOrDefault($v.CHROM, high(int))
      buf.add(BufferedRep(chromOrderIdx: chromOI, pos: v.POS, rec: bcf_dup(v.c)))

    vcf.close()

  # Sort by (chrom-order, pos).
  buf.sort(proc(a, b: BufferedRep): int =
    let c = cmp(a.chromOrderIdx, b.chromOrderIdx)
    if c != 0: c else: cmp(a.pos, b.pos))

  # Open output, write header, write records.
  let outPath = if isStdoutPath(cfg.outputPath): "/dev/stdout" else: cfg.outputPath
  let mode =
    if cfg.outputPath.endsWith(".bcf"):      "wb"
    elif cfg.outputPath.endsWith(".vcf.gz"): "wz"
    else:                                    "w"

  var outVcf: VCF
  if not open(outVcf, outPath, mode = mode):
    raise newException(IOError, "cannot open output: " & outPath)
  block:
    var dummy: VCF
    discard open(dummy, cfg.callers[0].path)
    outVcf.copy_header(dummy.header)
    dummy.close()
  bcf_hdr_destroy(outVcf.header.hdr)
  outVcf.header.hdr = finalHdr  # SHARED — clear before close to avoid double-free
  discard outVcf.write_header()

  for br in buf:
    discard bcf_write(vcfHtsFile(outVcf), finalHdr, br.rec)
    bcf_destroy(br.rec)

  outVcf.header.hdr = nil
  outVcf.close()
  logV("collapse: wrote " & $buf.len & " record(s) to " &
       (if isStdoutPath(cfg.outputPath): "stdout" else: cfg.outputPath))

# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

const CollapseVersion {.strdefine.} = "dev"

proc runCollapse*(cfg: CollapseConfig; cmdLine: string = "") =
  logV("matcha collapse: " & $cfg.callers.len & " caller(s)" &
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
  logV("resolving headers")
  let mh = resolveHeaders(cfg.callers)
  for w in mh.warnings: logWarn("collapse header: " & w)

  # First caller's chrom order anchors output ordering.
  var orderVcf: VCF
  if not open(orderVcf, cfg.callers[0].path):
    raise newException(IOError, "cannot open: " & cfg.callers[0].path)
  let chromOrder = captureChromOrder(orderVcf.header)
  orderVcf.close()

  # Phase 2: build shared output header, then integrated preproc+merge.
  var outSampleName: string
  let finalHdr = buildFinalHdr(cfg.callers, mh, cfg, CollapseVersion, cmdLine,
                               outSampleName)
  logV("integrated preproc+merge over " & $cfg.callers.len & " caller(s)")
  let msc = MergeStreamConfig(formatFields: cfg.formatFields,
                               nThreads:     cfg.nThreads,
                               tmpDir:       cfg.tmpDir)
  let im = integratedMerge(cfg.callers, mh, finalHdr, msc, chromOrder)

  # Build a PreprocOutput describing the merged slim BCFs for buildWorkQueue.
  let mergedPreproc = PreprocOutput(
    paths:          im.paths,
    populatedBins:  im.populatedBins,
    chromsBySvtype: im.chromsBySvtype,
    chromOrder:     chromOrder,
  )

  # Phase 3: Pass 1 — self-match over merged slim BCFs (emitSingletons=true).
  logV("self-matching merged slim BCFs")
  let matchCfg = MatchConfig(
    metric:         cfg.metric,
    threshold:      cfg.threshold,
    bndSlop:        cfg.bndSlop,
    nThreads:       cfg.nThreads,
    tmpDir:         cfg.tmpDir,
    selfMode:       true,
    emitSingletons: true,
  )
  let (jobs, fileList) = buildWorkQueue(mergedPreproc, mergedPreproc, matchCfg)
  logV("collapse self-match: " & $jobs.len & " job(s)")
  let allPairs = block:
    var r: seq[MatchPair]
    for jrs in runMatchPairJobsWithPool(jobs, matchCfg):
      for mp in jrs: r.add(mp)
    r
  logV("self-match: " & $allPairs.len & " pair(s)")

  # Phase 4: Build allOffsets from MatchPairs; build locByIdx for retrieval.
  let simMap = buildSimilarityMap(allPairs)
  var seenOffsets: HashSet[int32]
  var allOffsets: seq[int32]
  var locByIdx: Table[int32, tuple[chromIdx: int16, pos: int32, fileIdx: int16]]
  for p in allPairs:
    if p.srcIndexA notin seenOffsets:
      seenOffsets.incl(p.srcIndexA)
      allOffsets.add(p.srcIndexA)
      locByIdx[p.srcIndexA] = (p.chromIdx, p.posA, p.fileIdxA)
    if p.srcIndexB != NO_MATCH and p.srcIndexB notin seenOffsets:
      seenOffsets.incl(p.srcIndexB)
      allOffsets.add(p.srcIndexB)
      locByIdx[p.srcIndexB] = (p.chromIdx, p.posB, p.fileIdxB)
  logV("unique records: " & $allOffsets.len)

  # Phase 5: cluster.
  let clusters = clusterAll(allOffsets, simMap, cfg.linkage, cfg.threshold)

  # Build passQualMap via targeted CSI queries for all cluster members.
  var allMemberIdxs: HashSet[int32]
  for cl in clusters:
    for idx in cl: allMemberIdxs.incl(idx)
  var passQualMap: Table[int32, tuple[hasPASS: bool; qual: float32; callerIdx: int32]]
  var idxData, ciData: seq[int32]
  var membersByFile: Table[int16, seq[int32]]
  for idx in allMemberIdxs:
    let loc = locByIdx[idx]
    membersByFile.mgetOrPut(loc.fileIdx, @[]).add(idx)
  for fileIdx, members in membersByFile:
    var vcf: VCF
    if not open(vcf, fileList[fileIdx]):
      raise newException(IOError, "cannot open merged BCF: " & fileList[fileIdx])
    for idx in members:
      let loc = locByIdx[idx]
      let region = chromOrder[loc.chromIdx] & ":" & $loc.pos & "-" & $loc.pos
      for v in vcf.query(region):
        if readSrcIndex(v, idxData) != idx: continue
        let ci = if v.info().get("CALLER_IDX", ciData) == Status.OK and ciData.len > 0:
                   ciData[0] else: 0'i32
        passQualMap[idx] = (hasPASS: $v.FILTER == "PASS",
                            qual: v.QUAL.float32, callerIdx: ci)
        break
    vcf.close()

  # Select representatives.
  var finalClusters: seq[seq[int32]]
  for cl in clusters:
    let rep = selectRepresentative(cl, simMap, passQualMap, cfg.priority)
    var ordered = @[rep]
    for idx in cl:
      if idx != rep: ordered.add(idx)
    finalClusters.add(ordered)
  logV("clusters: " & $finalClusters.len)

  # Phase 6: write output.
  writeOutput(cfg, finalHdr, chromOrder, im.paths, finalClusters, passQualMap)

  # Clean up: merged slim BCFs + CSI indexes, and finalHdr.
  removeTempBcfs(im.paths)
  bcf_hdr_destroy(finalHdr)
