## collapse.nim — `matcha collapse` subcommand.
##
## Two-pass pipeline:
##   1. resolveHeaders — analyse N input headers, produce MergedHeader.
##   2. preprocessVcf (per caller, noIndex=true, keepPassQual=true) → caller slim BCFs.
##   3. mergeSortSlimBcfs — k-way merge per (svtype, bin) → merged slim BCFs.
##   4. runMatchJobsWithPool (self-mode) — Pass 1 matching over merged slim BCFs.
##   5. exploreMerged — Pass 2: enumerate allOffsets + passQualMap from merged BCFs.
##   6. clusterAll → selectRepresentative — cluster and pick representatives.
##   7. writeOutput — stream original caller files, buffer representative records,
##      sort by coordinate, write to output VCF/BCF.

import std/[algorithm, os, sets, strutils, tables]
import hts
import hts/private/hts_concat
import utils, preproc, match, matchcore, mergecore, log

from hts/private/hts_concat import libname

proc bcf_translate*(dst_hdr, src_hdr: ptr bcf_hdr_t; line: ptr bcf1_t): cint
  {.cdecl, importc: "bcf_translate", dynlib: libname.}

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
    infoFields*:   seq[string]   ## --info filter; empty = keep all
    outputPath*:   string
    nThreads*:     int
    tmpDir*:       string
    callers*:      seq[CallerInput]

# ---------------------------------------------------------------------------
# Header traversal helpers (reused from mergecore context)
# ---------------------------------------------------------------------------

proc collectFilterLines(h: ptr bcf_hdr_t): seq[string] =
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_FLT.cint: continue
    let keys = cast[ptr UncheckedArray[cstring]](hr.keys)
    let vals = cast[ptr UncheckedArray[cstring]](hr.vals)
    var line = "##FILTER=<"
    for j in 0 ..< hr.nkeys.int:
      if j > 0: line &= ","
      line &= $keys[j] & "=" & $vals[j]
    line &= ">"
    result.add(line)

proc infoFieldDefs(h: ptr bcf_hdr_t): seq[(string, string, string, string)] =
  ## Return (id, number, type, desc) for each INFO field in the header.
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

proc infoFieldNames(h: ptr bcf_hdr_t): seq[string] =
  for (id, _, _, _) in infoFieldDefs(h): result.add(id)

proc fmtFieldNames(h: ptr bcf_hdr_t): seq[string] =
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_FMT.cint: continue
    let keys = cast[ptr UncheckedArray[cstring]](hr.keys)
    let vals = cast[ptr UncheckedArray[cstring]](hr.vals)
    for j in 0 ..< hr.nkeys.int:
      if $keys[j] == "ID": result.add($vals[j]); break

# ---------------------------------------------------------------------------
# Pass 2 — enumerate allOffsets + PASS/QUAL from merged slim BCFs
# ---------------------------------------------------------------------------

proc exploreMerged(mergedPaths: Table[SvtypeBin, string]):
    tuple[allOffsets: seq[int64];
          passQualMap: Table[int64, tuple[hasPASS: bool; qual: float32]]] =
  var seen: HashSet[int64]
  var boffScratch: seq[int32]
  for path in mergedPaths.values:
    var vcf: VCF
    if not open(vcf, path):
      raise newException(IOError, "cannot open merged slim BCF: " & path)
    for v in vcf:
      let off = readBoff(v, boffScratch)
      if off in seen: continue
      seen.incl(off)
      result.allOffsets.add(off)
      let hasPASS = $v.FILTER == "PASS"
      let qual = v.QUAL.float32
      result.passQualMap[off] = (hasPASS: hasPASS, qual: qual)
    vcf.close()

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

type BufferedRep = object
  chromOrderIdx: int
  pos:           int64
  rec:           ptr bcf1_t  ## owned; caller must bcf_destroy

proc buildOutputHdr(cfg: CollapseConfig; mh: MergedHeader;
                    version, cmdLine: string): ptr bcf_hdr_t =
  ## Build output header: first caller contigs + all callers' contigs,
  ## merged INFO/FORMAT (minus MATCHA_CALLER_IDX, filtered by --info/--format),
  ## SOURCE/SOURCELIST/N_SOURCE/N_MERGED, provenance comments.
  var firstVcf: VCF
  if not open(firstVcf, cfg.callers[0].path):
    raise newException(IOError, "cannot open: " & cfg.callers[0].path)
  result = bcf_hdr_dup(firstVcf.header.hdr)
  firstVcf.close()

  bcf_hdr_remove(result, BCF_HEADER_TYPE.BCF_HL_INFO.cint, nil)
  bcf_hdr_remove(result, BCF_HEADER_TYPE.BCF_HL_FMT.cint, nil)

  # Add merged defs — skip internal field, skip INFO fields not in --info filter.
  let infoFilter = toHashSet(cfg.infoFields)
  proc keepInfoField(name: string): bool =
    if name in ["MATCHA_CALLER_IDX", "MATCHA_BOFF"]: return false
    if infoFilter.len == 0: return true
    for tok in cfg.infoFields:
      if name == tok or name.startsWith(tok & "_"): return true
    false

  for line in mh.headerLines:
    # Extract the ID from the header line to check the filter.
    let idStart = line.find("ID=")
    if idStart < 0:
      discard bcf_hdr_append(result, line.cstring)
      continue
    let idEnd = line.find(',', idStart + 3)
    let fieldId = if idEnd > 0: line[idStart + 3 ..< idEnd] else: line[idStart + 3 .. ^2]
    let isFmt = line.startsWith("##FORMAT")
    let fmtKeep = toHashSet(cfg.formatFields)
    if isFmt:
      if fieldId in fmtKeep:
        discard bcf_hdr_append(result, line.cstring)
    else:
      if keepInfoField(fieldId):
        discard bcf_hdr_append(result, line.cstring)

  # FILTER defs: merge from all callers.
  var seenFlt: HashSet[string]
  for caller in cfg.callers:
    var vcf: VCF
    if not open(vcf, caller.path): continue
    for line in collectFilterLines(vcf.header.hdr):
      if line notin seenFlt:
        seenFlt.incl(line)
        discard bcf_hdr_append(result, line.cstring)
    vcf.close()

  # Contigs from callers 1+ (first caller's contigs are already in the dup'd header).
  if cfg.callers.len > 1:
    var seenCtg: HashSet[string]
    var firstVcf2: VCF
    if open(firstVcf2, cfg.callers[0].path):
      var n: cint = 0
      let names = bcf_hdr_seqnames(firstVcf2.header.hdr, n.addr)
      if names != nil:
        for i in 0 ..< n.int: seenCtg.incl($names[i])
        free(names)
      firstVcf2.close()
    for ci in 1 ..< cfg.callers.len:
      var vcf: VCF
      if not open(vcf, cfg.callers[ci].path): continue
      var n: cint = 0
      let names = bcf_hdr_seqnames(vcf.header.hdr, n.addr)
      if names != nil:
        for i in 0 ..< n.int:
          let c = $names[i]
          if c notin seenCtg:
            seenCtg.incl(c)
            discard bcf_hdr_append(result, ("##contig=<ID=" & c & ">").cstring)
        free(names)
      vcf.close()

  # Provenance INFO fields.
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
  discard bcf_hdr_sync(result)

proc isStdoutPath(p: string): bool = p == "" or p == "-" or p == "/dev/stdout"

type ClusterProv = object
  source:     string
  sourceList: seq[string]
  nSource:    int
  nMerged:    int

proc writeOutput(cfg: CollapseConfig;
                  outHdr: ptr bcf_hdr_t;
                  chromOrder: seq[string];
                  finalClusters: seq[seq[int64]];
                  mh: MergedHeader) =

  # Build provenance map: repCompositeOff → ClusterProv.
  # callerIdx is decoded from the high 16 bits of each compositeOff.
  var repProv: Table[int64, ClusterProv]
  for cl in finalClusters:
    if cl.len == 0: continue
    let repOff = cl[0]
    var callerIdxSeen: seq[int]
    for off in cl:
      let ci = int(off shr 48)
      if ci notin callerIdxSeen: callerIdxSeen.add(ci)
    callerIdxSeen.sort()
    var sourceList: seq[string]
    for ci in callerIdxSeen: sourceList.add(cfg.callers[ci].name)
    repProv[repOff] = ClusterProv(
      source:     cfg.callers[int(repOff shr 48)].name,
      sourceList: sourceList,
      nSource:    sourceList.len,
      nMerged:    cl.len,
    )

  # Build chrom→order index.
  var chromIdx: Table[string, int]
  for i, c in chromOrder.pairs: chromIdx[c] = i

  # INFO field filter predicate (closure).
  let infoFilterSet = cfg.infoFields
  proc keepInfo(name: string): bool =
    if name in ["MATCHA_CALLER_IDX", "MATCHA_BOFF"]: return false
    if infoFilterSet.len == 0: return true
    for tok in infoFilterSet:
      if name == tok or name.startsWith(tok & "_"): return true
    false

  let fmtKeepSet = toHashSet(cfg.formatFields)

  # Group representatives by callerIdx → set of origOffs.
  # compositeOff = (callerIdx shl 48) or origOff
  var repByCallerOrig: Table[int, HashSet[int64]]
  for repOff in repProv.keys:
    let ci = int(repOff shr 48)
    let origOff = repOff and 0x0000_FFFF_FFFF_FFFF'i64
    repByCallerOrig.mgetOrPut(ci, initHashSet[int64]()).incl(origOff)

  # Stream each original caller file, collect representative records.
  var buf: seq[BufferedRep]
  for ci, caller in cfg.callers.pairs:
    if ci notin repByCallerOrig: continue
    let origOffsets = repByCallerOrig[ci]

    var vcf: VCF
    if not open(vcf, caller.path):
      raise newException(IOError, "cannot open: " & caller.path)

    # Add renamed INFO defs to source header (needed for in-place rename + translate).
    for origName, res in mh.infoRes.pairs:
      if res.kind == fcIncompatibleInfo and ci in res.renames:
        let newName = res.renames[ci]
        for (id, num, typ, desc) in infoFieldDefs(vcf.header.hdr):
          if id == origName:
            discard vcf.header.add_info(newName, num, typ, desc); break

    # Add provenance fields to source header (needed for v.info.set).
    discard vcf.header.add_info("SOURCE",     "1", "String",  "Representative caller name")
    discard vcf.header.add_info("SOURCELIST", ".", "String",  "All caller names in cluster")
    discard vcf.header.add_info("N_SOURCE",   "1", "Integer", "Distinct input callsets")
    discard vcf.header.add_info("N_MERGED",   "1", "Integer", "Total records merged")

    var nextOff = int64(bgzf_tell(bgzfHandle(vcf)))
    for v in vcf:
      let recOff = nextOff
      nextOff = int64(bgzf_tell(bgzfHandle(vcf)))

      if recOff notin origOffsets: continue

      # Recompose the offset and look up provenance in O(1).
      let compositeOff = (ci.int64 shl 48) or recOff
      if compositeOff notin repProv: continue
      let prov = repProv[compositeOff]

      # Apply INFO renames.
      for origName, res in mh.infoRes.pairs:
        if res.kind != fcIncompatibleInfo: continue
        let newName = res.renames.getOrDefault(ci, "")
        if newName.len == 0: continue
        var intData: seq[int32]
        var fltData: seq[float32]
        var strData: string
        if v.info().get(origName, intData) == Status.OK:
          discard v.info.set(newName, intData)
        elif v.info().get(origName, fltData) == Status.OK:
          discard v.info.set(newName, fltData)
        elif v.info().get(origName, strData) == Status.OK:
          discard v.info.set(newName, strData)
        discard v.info.delete(origName)

      # Drop INFO fields not passing the filter.
      var toDeleteInfo: seq[string]
      for fld in v.info.fields:
        if not keepInfo(fld.name): toDeleteInfo.add(fld.name)
      for name in toDeleteInfo:
        discard v.info.delete(name)

      # Drop FORMAT fields not in fmtKeep.
      var fmtNames: seq[string]
      for fld in v.format.fields: fmtNames.add(fld.name)
      for name in fmtNames:
        if name notin fmtKeepSet:
          try: discard v.format.delete(name)
          except KeyError: discard

      # Add provenance INFO.
      var srcStr = prov.source
      discard v.info.set("SOURCE", srcStr)
      var nSrc = prov.nSource.int32
      discard v.info.set("N_SOURCE", nSrc)
      var nMrg = prov.nMerged.int32
      discard v.info.set("N_MERGED", nMrg)
      var slStr = prov.sourceList.join(",")
      discard v.info.set("SOURCELIST", slStr)

      # Translate field IDs into output header space, then copy.
      discard bcf_translate(outHdr, vcf.header.hdr, v.c)
      let chromOI = chromIdx.getOrDefault($v.CHROM, high(int))
      buf.add(BufferedRep(chromOrderIdx: chromOI, pos: v.POS, rec: bcf_dup(v.c)))

    vcf.close()

  # Sort by coordinate.
  buf.sort(proc(a, b: BufferedRep): int =
    let c = cmp(a.chromOrderIdx, b.chromOrderIdx)
    if c != 0: c else: cmp(a.pos, b.pos))

  # Open output and write.
  let outPath = if isStdoutPath(cfg.outputPath): "/dev/stdout" else: cfg.outputPath
  let mode =
    if cfg.outputPath.endsWith(".bcf"):    "wb"
    elif cfg.outputPath.endsWith(".vcf.gz"): "wz"
    else:                                   "w"

  var outVcf: VCF
  if not open(outVcf, outPath, mode = mode):
    raise newException(IOError, "cannot open output: " & outPath)
  # open() in write mode leaves outVcf.header nil; copy_header initialises it.
  block:
    var srcVcf: VCF
    discard open(srcVcf, cfg.callers[0].path)
    outVcf.copy_header(srcVcf.header)
    srcVcf.close()
  bcf_hdr_destroy(outVcf.header.hdr)
  outVcf.header.hdr = outHdr
  discard outVcf.write_header()

  for br in buf:
    discard bcf_write(vcfHtsFile(outVcf), outHdr, br.rec)
    bcf_destroy(br.rec)

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
      stderr.writeLine "error: input file not found: " & caller.path
      quit(1)
  if cfg.tmpDir == "":
    stderr.writeLine "error: tmpDir must be set"
    quit(1)

  # Phase 1: resolve output header from all N input headers.
  logV("resolving headers")
  let mh = resolveHeaders(cfg.callers)
  for w in mh.warnings: logWarn("collapse header: " & w)

  # Use the first caller's chrom order throughout (preproc/match/output).
  var orderVcf: VCF
  if not open(orderVcf, cfg.callers[0].path):
    raise newException(IOError, "cannot open: " & cfg.callers[0].path)
  let chromOrder = captureChromOrder(orderVcf.header)
  orderVcf.close()

  let ioThreads = if cfg.nThreads >= 2: 2 else: 0

  # Phase 2: preprocess each caller into per-(svtype, bin) slim BCFs.
  # noIndex=true: CSI index built later on the merged output only.
  # keepPassQual=true: preserve QUAL for representative selection in Pass 2.
  logV("preprocessing " & $cfg.callers.len & " caller(s)")
  var callerOutputs = newSeq[PreprocOutput](cfg.callers.len)
  for ci, caller in cfg.callers.pairs:
    callerOutputs[ci] = preprocessVcf(
      caller.path, cfg.tmpDir, caller.name,
      ioThreads = ioThreads, noIndex = true, keepPassQual = true)

  # Phase 3: merge-sort per-(svtype, bin) slim BCFs → one merged slim BCF per key.
  logV("merge-sorting slim BCFs")
  var mergedPaths: Table[SvtypeBin, string]
  var mergedPop:   Table[SvType, set[uint8]]
  var mergedChroms: Table[SvType, HashSet[string]]

  var allKeys: HashSet[SvtypeBin]
  for co in callerOutputs:
    for key in co.paths.keys: allKeys.incl(key)
    for svt, bins in co.populatedBins.pairs:
      if svt notin mergedPop: mergedPop[svt] = {}
      mergedPop[svt] = mergedPop[svt] + bins
    for svt, chroms in co.chromsBySvtype.pairs:
      if svt notin mergedChroms: mergedChroms[svt] = initHashSet[string]()
      for c in chroms: mergedChroms[svt].incl(c)

  for key in allKeys:
    var inputs: seq[string]
    var callerIdxs: seq[int]
    for ci, co in callerOutputs.pairs:
      if key in co.paths:
        inputs.add(co.paths[key])
        callerIdxs.add(ci)
    let (svt, bin) = key
    let outPath = cfg.tmpDir / "matcha_" & $getCurrentProcessId() &
                  "_M_" & $svt & "_b" & $bin & ".bcf"
    mergeSortSlimBcfs(inputs, callerIdxs, outPath, chromOrder)
    mergedPaths[key] = outPath

  # Build a PreprocOutput describing the merged slim BCFs for buildWorkQueue.
  var mergedPreproc = PreprocOutput(
    paths:          mergedPaths,
    populatedBins:  mergedPop,
    chromsBySvtype: mergedChroms,
    chromOrder:     chromOrder,
  )

  # Phase 4: Pass 1 — self-match over merged slim BCFs.
  logV("self-matching merged slim BCFs")
  let matchCfg = MatchConfig(
    metric:    cfg.metric,
    threshold: cfg.threshold,
    bndSlop:   cfg.bndSlop,
    nThreads:  cfg.nThreads,
    tmpDir:    cfg.tmpDir,
    selfMode:  true,
  )
  let jobs = buildWorkQueue(mergedPreproc, mergedPreproc, matchCfg)
  logV("collapse self-match: " & $jobs.len & " job(s)")
  let allPairs = block:
    var r: seq[MatchPair]
    for jrs in runMatchPairJobsWithPool(jobs, matchCfg):
      for mp in jrs: r.add(mp)
    r
  logV("self-match: " & $allPairs.len & " pair(s)")

  # Phase 5: Pass 2 — enumerate allOffsets + PASS/QUAL from merged slim BCFs.
  let simMap = buildSimilarityMap(allPairs)
  let (allOffsets, passQualMap) = exploreMerged(mergedPaths)
  logV("offsets seen: " & $allOffsets.len)

  # Phase 6: cluster and select representatives.
  let clusters = clusterAll(allOffsets, simMap, cfg.linkage, cfg.threshold)
  var finalClusters: seq[seq[int64]]
  for cl in clusters:
    let rep = selectRepresentative(cl, simMap, passQualMap, cfg.priority)
    var ordered = @[rep]
    for off in cl:
      if off != rep: ordered.add(off)
    finalClusters.add(ordered)
  logV("clusters: " & $finalClusters.len)

  # Phase 7: write output.
  let outHdr = buildOutputHdr(cfg, mh, CollapseVersion, cmdLine)
  writeOutput(cfg, outHdr, chromOrder, finalClusters, mh)

  # Clean up: merged slim BCFs + per-caller slim BCFs.
  for path in mergedPaths.values:
    if fileExists(path):          removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")
  for co in callerOutputs:
    for path in co.paths.values:
      if fileExists(path): removeFile(path)
