## collapse.nim — `matcha collapse` subcommand.
##
## Single-pass pipeline (fused preproc + merge via synced_bcf_reader):
##   1. resolveHeaders — analyse N input headers, produce MergedHeader.
##   2. integratedMerge — stream all N caller VCFs in lockstep via one
##      `bcf_srs_t`, normalize + filter + write per-(svtype, bin) merged BCFs
##      in a single pass. Uses a shared htsThreadPool for parallel BGZF
##      decompression across readers + writers.
##   3. runMatchJobsWithPool (self-mode) — Pass 1 matching over merged slim BCFs.
##   4. exploreMerged — Pass 2: enumerate allOffsets + passQualMap.
##   5. clusterAll → selectRepresentative — cluster and pick representatives.
##   6. writeOutput — stream merged BCFs, buffer representative records,
##      apply output-time INFO filter, sort by coordinate, write final VCF/BCF.
##
## MATCHA_BOFF in collapse is an opaque identity token: `(callerIdx << 48) |
## monotonic_counter`. It is no longer a real BGZF offset. matchcore's
## `selectRepresentative` still works (high 16 bits = caller index).

import std/[algorithm, os, sets, strutils, tables]
import hts
import hts/private/hts_concat
import utils, preproc, match, matchcore, mergecore, log, synced_bcf_reader

from hts/private/hts_concat import libname

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

proc fmtFieldDefs(h: ptr bcf_hdr_t): seq[(string, string, string, string)] =
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_FMT.cint: continue
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
# integratedMerge — fused preproc + merge in one synced_bcf_reader pass
# ---------------------------------------------------------------------------

const AlwaysKeepInMerged = ["SVTYPE", "SVLEN", "END", "CHR2", "POS2",
                            "SOURCE", "SOURCELIST", "N_SOURCE", "N_MERGED",
                            "MATCHA_BOFF"]

proc keepInfoForMerged(name: string; infoFilter: seq[string]): bool =
  ## Keep set for INFO fields in the merged slim BCFs.
  ## - Always keeps matcha-internal + matchcore-required fields.
  ## - With infoFilter empty: keep all non-internal fields.
  ## - With infoFilter set: keep listed fields (or their *_<caller> renames).
  if name == "MATCHA_CALLER_IDX": return false  # legacy internal
  for n in AlwaysKeepInMerged:
    if name == n: return true
  if infoFilter.len == 0: return true
  for tok in infoFilter:
    if name == tok or name.startsWith(tok & "_"): return true
  false

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
      if keepInfoForMerged(fieldId, cfg.infoFields):
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

  # Provenance / matcha-internal INFO defs.
  discard bcf_hdr_append(result,
    "##INFO=<ID=MATCHA_BOFF,Number=2,Type=Integer,Description=\"matcha-internal: composite identity token (callerIdx<<48 | counter)\">".cstring)
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

proc augmentSrcHdrForRenames(srcHdr: ptr bcf_hdr_t; ci: int;
                              mh: MergedHeader) =
  ## For each renamed INFO/FORMAT field this caller writes under, append the
  ## renamed def to srcHdr so subsequent bcf_update_info/_format calls have
  ## the field's Number/Type available. Must be done before any rename write.
  for origName, res in mh.infoRes.pairs:
    if res.kind != fcIncompatibleInfo: continue
    if ci notin res.renames: continue
    let newName = res.renames[ci]
    for (id, num, typ, desc) in infoFieldDefs(srcHdr):
      if id == origName:
        discard bcf_hdr_append(srcHdr,
          ("##INFO=<ID=" & newName & ",Number=" & num & ",Type=" & typ &
           ",Description=\"" & desc & "\">").cstring)
        break
  for origName, res in mh.fmtRes.pairs:
    if res.kind != fcIncompatibleFmt: continue
    if ci notin res.renames: continue
    let newName = res.renames[ci]
    for (id, num, typ, desc) in fmtFieldDefs(srcHdr):
      if id == origName:
        discard bcf_hdr_append(srcHdr,
          ("##FORMAT=<ID=" & newName & ",Number=" & num & ",Type=" & typ &
           ",Description=\"" & desc & "\">").cstring)
        break
  # Also ensure SOURCE / MATCHA_BOFF / SV defs exist for the authoritative
  # writes we issue against srcHdr in the per-record loop.
  proc ensureInfo(h: ptr bcf_hdr_t; name, num, typ, desc: string) =
    for (id, _, _, _) in infoFieldDefs(h):
      if id == name: return
    discard bcf_hdr_append(h,
      ("##INFO=<ID=" & name & ",Number=" & num & ",Type=" & typ &
       ",Description=\"" & desc & "\">").cstring)
  ensureInfo(srcHdr, "MATCHA_BOFF", "2", "Integer", "matcha-internal identity token")
  ensureInfo(srcHdr, "SOURCE",      "1", "String",  "Representative caller name")
  ensureInfo(srcHdr, "SVTYPE",      "1", "String",  "Type of structural variant")
  ensureInfo(srcHdr, "SVLEN",       "1", "Integer", "Length of the SV")
  ensureInfo(srcHdr, "END",         "1", "Integer", "End position of the SV")
  ensureInfo(srcHdr, "CHR2",        "1", "String",  "Chromosome of mate breakend")
  ensureInfo(srcHdr, "POS2",        "1", "Integer", "Position of mate breakend")
  discard bcf_hdr_sync(srcHdr)

proc btToHt(bt: cint): cint =
  if bt == BCF_BT_FLOAT: BCF_HT_REAL.cint
  elif bt == BCF_BT_CHAR: BCF_HT_STR.cint
  else: BCF_HT_INT.cint

proc applyInfoRename(srcHdr: ptr bcf_hdr_t; rec: ptr bcf1_t;
                     origName, newName: string;
                     buf: var pointer; bufN: var cint) =
  let info = bcf_get_info(srcHdr, rec, origName.cstring)
  if info == nil: return
  let htType = btToHt(info.`type`)
  let n = bcf_get_info_values(srcHdr, rec, origName.cstring,
                              buf.addr, bufN.addr, htType)
  if n > 0:
    discard bcf_update_info(srcHdr, rec, newName.cstring, buf, n.cint, htType)
  discard bcf_update_info(srcHdr, rec, origName.cstring, nil, 0.cint, htType)

proc applyFmtRename(srcHdr: ptr bcf_hdr_t; rec: ptr bcf1_t;
                    origName, newName: string;
                    buf: var pointer; bufN: var cint) =
  let fmt = bcf_get_fmt(srcHdr, rec, origName.cstring)
  if fmt == nil: return
  let htType = btToHt(fmt.`type`)
  let n = bcf_get_format_values(srcHdr, rec, origName.cstring,
                                buf.addr, bufN.addr, htType)
  if n > 0:
    discard bcf_update_format(srcHdr, rec, newName.cstring, buf, n.cint, htType)
  discard bcf_update_format(srcHdr, rec, origName.cstring, nil, 0.cint, htType)

type IntegratedMergeResult* = object
  paths*:          Table[SvtypeBin, string]
  populatedBins*:  Table[SvType, set[uint8]]
  chromsBySvtype*: Table[SvType, HashSet[string]]
  finalHdr*:       ptr bcf_hdr_t

proc integratedMerge*(callers: seq[CallerInput]; mh: MergedHeader;
                     chromOrder: seq[string]; cfg: CollapseConfig;
                     version, cmdLine: string): IntegratedMergeResult =
  ## Stream all N caller VCFs via one synced_bcf_reader, normalize each
  ## record, filter INFO/FORMAT to the user-selected fields, write per-
  ## (svtype, bin) merged slim BCFs. Returns the merged paths + metadata
  ## + the shared finalHdr used for both merged BCFs and final output.

  # 1. Build the shared output header.
  var outSampleName: string
  result.finalHdr = buildFinalHdr(callers, mh, cfg, version, cmdLine, outSampleName)

  # 2. Init synced_bcf_reader + thread pool.
  let sr = bcf_sr_init()
  if sr == nil:
    raise newException(IOError, "bcf_sr_init failed")
  discard bcf_sr_set_opt(sr, BCF_SR_ALLOW_NO_IDX)

  let nThr = max(1, cfg.nThreads)
  var tpool = htsThreadPool(pool: nil, qsize: 0)
  if nThr >= 2:
    tpool.pool = hts_tpool_init(nThr.cint)

  for caller in callers:
    if bcf_sr_add_reader(sr, caller.path.cstring) == 0:
      let msg = $bcf_sr_strerror(srs_errnum(sr))
      bcf_sr_destroy(sr)
      if tpool.pool != nil: hts_tpool_destroy(tpool.pool)
      raise newException(IOError, "cannot open: " & caller.path & " (" & msg & ")")

  let nreaders = srs_nreaders(sr).int

  # Attach thread pool to each reader's underlying file.
  if tpool.pool != nil:
    for i in 0 ..< nreaders:
      discard hts_set_opt(srs_get_file(sr, i.cint),
                          srs_hts_opt_thread_pool(), tpool.addr)

  # Subset samples on each reader's header (one sample per caller for
  # collapse, or none).
  for ci in 0 ..< nreaders:
    let srcHdr = srs_get_header(sr, ci.cint)
    if cfg.formatFields.len == 0:
      discard bcf_hdr_set_samples(srcHdr, cstring(nil), 0.cint)
    else:
      let nsamp = bcf_hdr_nsamples(srcHdr).int
      if nsamp > 0:
        let firstSample = $cast[cstringArray](srcHdr.samples)[0]
        discard bcf_hdr_set_samples(srcHdr, firstSample.cstring, 0.cint)
      else:
        discard bcf_hdr_set_samples(srcHdr, cstring(nil), 0.cint)
    augmentSrcHdrForRenames(srcHdr, ci, mh)

  # Per-caller WarnState + record counter.
  var wsList: seq[WarnState]
  var recCounter = newSeq[int64](nreaders)
  for i in 0 ..< nreaders:
    wsList.add(initWarnState(callers[i].name))

  # Keep sets for filtering, derived from finalHdr (so post-filter records
  # contain exactly the fields finalHdr defines → bcf_translate is clean).
  var infoKeepSet, fmtKeepSet: HashSet[string]
  for (id, _, _, _) in infoFieldDefs(result.finalHdr): infoKeepSet.incl(id)
  for (id, _, _, _) in fmtFieldDefs(result.finalHdr):  fmtKeepSet.incl(id)

  # Reusable Variant view over the synced reader's records; one per loop,
  # re-pointed at each record. See newVariantView for lifetime notes.
  let view = newVariantView()
  var svtypeBuf: string
  var endBuf, svlenBuf: seq[int32]
  var renameBuf: pointer = nil
  var renameN:   cint    = 0

  # Lazy-opened writers + CSI indexes per (svtype, bin).
  var writers:    Table[SvtypeBin, VCF]
  var indexes:    Table[SvtypeBin, ptr hts_idx_t]
  var writerHdrs: Table[SvtypeBin, ptr bcf_hdr_t]

  try:
    # 3. Stream records.
    while bcf_sr_next_line(sr) > 0:
      for ci in 0 ..< nreaders:
        if srs_has_line(sr, ci.cint) == 0: continue
        let srcHdr = srs_get_header(sr, ci.cint)
        let rawRec = srs_get_line(sr, ci.cint)

        let counter = recCounter[ci]
        inc recCounter[ci]

        let rec = bcf_dup(rawRec)

        let nr = normalizeRecord(srcHdr, rec, counter.int + 1,
                                 wsList[ci], view, svtypeBuf, endBuf, svlenBuf,
                                 callers[ci].path)
        if not nr.ok:
          bcf_destroy(rec)
          continue

        # 3a. Apply INFO renames.
        for origName, res in mh.infoRes.pairs:
          if res.kind != fcIncompatibleInfo: continue
          if ci notin res.renames: continue
          applyInfoRename(srcHdr, rec, origName, res.renames[ci],
                          renameBuf, renameN)

        # 3b. Apply FORMAT renames.
        for origName, res in mh.fmtRes.pairs:
          if res.kind != fcIncompatibleFmt: continue
          if ci notin res.renames: continue
          applyFmtRename(srcHdr, rec, origName, res.renames[ci],
                         renameBuf, renameN)

        # 3c. Filter INFO fields: drop any not in infoKeepSet.
        discard bcf_unpack(rec, BCF_UN_INFO.cint)
        block:
          let nInfo = rec.n_info.int
          if nInfo > 0:
            var toDel: seq[string]
            let infoArr = cast[ptr UncheckedArray[bcf_info_t]](rec.d.info)
            let idPairs = cast[ptr UncheckedArray[bcf_idpair_t]](srcHdr.id[BCF_DT_ID])
            for i in 0 ..< nInfo:
              let key = infoArr[i].key
              if key < 0: continue
              let nm = $idPairs[key].key
              if nm notin infoKeepSet: toDel.add(nm)
            for nm in toDel:
              discard bcf_update_info(srcHdr, rec, nm.cstring,
                                      nil, 0.cint, BCF_HT_INT.cint)

        # 3d. Filter FORMAT fields.
        discard bcf_unpack(rec, BCF_UN_FMT.cint)
        block:
          let nFmt = rec.n_fmt.int
          if nFmt > 0:
            var toDel: seq[string]
            let fmtArr = cast[ptr UncheckedArray[bcf_fmt_t]](rec.d.fmt)
            let idPairs = cast[ptr UncheckedArray[bcf_idpair_t]](srcHdr.id[BCF_DT_ID])
            for i in 0 ..< nFmt:
              let id = fmtArr[i].id
              if id < 0: continue
              let nm = $idPairs[id].key
              if nm notin fmtKeepSet: toDel.add(nm)
            for nm in toDel:
              discard bcf_update_format(srcHdr, rec, nm.cstring,
                                        nil, 0.cint, BCF_HT_INT.cint)

        # 3e. Authoritative writes for END / CHR2 / POS2 / SVTYPE / SVLEN.
        block:
          var svtStr = $nr.svt
          discard bcf_update_info(srcHdr, rec, "SVTYPE".cstring,
                                  svtStr[0].addr, svtStr.len.cint, BCF_HT_STR.cint)
          var svlenVal = nr.svlen.int32
          discard bcf_update_info(srcHdr, rec, "SVLEN".cstring,
                                  svlenVal.addr, 1.cint, BCF_HT_INT.cint)
          if nr.svt == svBND:
            var chr2Str = nr.bndChr2
            if chr2Str.len > 0:
              discard bcf_update_info(srcHdr, rec, "CHR2".cstring,
                                      chr2Str[0].addr, chr2Str.len.cint, BCF_HT_STR.cint)
            var pos2Val = nr.bndPos2.int32
            discard bcf_update_info(srcHdr, rec, "POS2".cstring,
                                    pos2Val.addr, 1.cint, BCF_HT_INT.cint)
          else:
            var endVal = nr.endPos.int32
            discard bcf_update_info(srcHdr, rec, "END".cstring,
                                    endVal.addr, 1.cint, BCF_HT_INT.cint)

        # 3f. Set SOURCE = caller name.
        block:
          var srcStr = callers[ci].name
          discard bcf_update_info(srcHdr, rec, "SOURCE".cstring,
                                  srcStr[0].addr, srcStr.len.cint, BCF_HT_STR.cint)

        # 3g. Set MATCHA_BOFF = (callerIdx << 48) | counter.
        block:
          let compositeOff = (ci.int64 shl 48) or counter
          var boffPair = encodeBoff(compositeOff)
          discard bcf_update_info(srcHdr, rec, "MATCHA_BOFF".cstring,
                                  boffPair[0].addr, 2.cint, BCF_HT_INT.cint)

        # 3h. REF/ALT trim to keep records small (matchcore reads neither).
        # NOTE: bcf_update_alleles_str expects comma-separated alleles.
        discard bcf_update_alleles_str(srcHdr, rec, "N,.".cstring)

        # 3i. Track metadata + lazy-open writer.
        if nr.svt notin result.populatedBins:
          result.populatedBins[nr.svt] = {}
        result.populatedBins[nr.svt].incl(uint8(nr.binIdx))
        if nr.svt notin result.chromsBySvtype:
          result.chromsBySvtype[nr.svt] = initHashSet[string]()
        result.chromsBySvtype[nr.svt].incl(getChromName(srcHdr, rec.rid))

        let key: SvtypeBin = (nr.svt, nr.binIdx)
        if key notin writers:
          let outPath = cfg.tmpDir / "matcha_" & $getCurrentProcessId() &
                        "_M_" & $nr.svt & "_b" & $nr.binIdx & ".bcf"
          var wtr: VCF
          if not open(wtr, outPath, mode = "wb"):
            raise newException(IOError, "cannot create merged BCF: " & outPath)
          # Init wtr.header via a dummy source, then replace with dup of finalHdr
          # (each writer owns its dup; finalHdr remains valid for writeOutput).
          var dummy: VCF
          discard open(dummy, callers[0].path)
          wtr.copy_header(dummy.header)
          dummy.close()
          bcf_hdr_destroy(wtr.header.hdr)
          let wHdr = bcf_hdr_dup(result.finalHdr)
          wtr.header.hdr = wHdr
          writerHdrs[key] = wHdr
          if tpool.pool != nil:
            discard hts_set_opt(vcfHtsFile(wtr),
                                srs_hts_opt_thread_pool(), tpool.addr)
          discard wtr.write_header()
          writers[key] = wtr
          result.paths[key] = outPath
          let headerOff = uint64(bgzf_tell(bgzfHandle(wtr)))
          let idx = hts_idx_init(0.cint, HTS_FMT_CSI.cint, headerOff,
                                 14.cint, 5.cint)
          if idx == nil:
            raise newException(IOError, "cannot create CSI index: " & outPath)
          indexes[key] = idx

        # 3j. Translate field IDs to finalHdr space, write record, push to index.
        discard bcf_translate(writerHdrs[key], srcHdr, rec)
        discard bcf_write(vcfHtsFile(writers[key]), writerHdrs[key], rec)
        let woff = uint64(bgzf_tell(bgzfHandle(writers[key])))
        discard hts_idx_push(indexes[key], rec.rid, int64(rec.pos),
                             nr.endPos, woff, 1.cint)

        bcf_destroy(rec)

  finally:
    # Teardown order: destroy synced reader first (it owns reader file handles),
    # then close writers (each holds its own dup of finalHdr), then destroy
    # the shared thread pool last (after every htsFile* that referenced it).
    bcf_sr_destroy(sr)
    for key, wtr in writers.mpairs:
      let finalOff = uint64(bgzf_tell(bgzfHandle(wtr)))
      hts_idx_finish(indexes[key], finalOff)
      wtr.close()  # bcf_hdr_destroy on wtr.header.hdr (the dup) — finalHdr is safe
      let path = result.paths[key]
      hts_idx_save(indexes[key], path.cstring, HTS_FMT_CSI.cint)
      hts_idx_destroy(indexes[key])
    if tpool.pool != nil: hts_tpool_destroy(tpool.pool)

    if renameBuf != nil: c_free(renameBuf)

  # Emit per-caller summaries.
  for ws in wsList: emitSummary(ws)

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

type
  BufferedRep = object
    chromOrderIdx: int
    pos:           int64
    rec:           ptr bcf1_t  ## owned; caller must bcf_destroy

  ClusterProv = object
    source:     string
    sourceList: seq[string]
    nSource:    int
    nMerged:    int

proc isStdoutPath(p: string): bool = p == "" or p == "-" or p == "/dev/stdout"

proc writeOutput(cfg: CollapseConfig;
                 finalHdr: ptr bcf_hdr_t;
                 chromOrder: seq[string];
                 mergedPaths: Table[SvtypeBin, string];
                 finalClusters: seq[seq[int64]]) =
  ## Stream merged BCFs, pick records whose MATCHA_BOFF marks them as cluster
  ## representatives, set provenance fields, apply output-time INFO filter,
  ## sort by coordinate, write final VCF/BCF. Records are already in finalHdr
  ## field-ID space (merged BCFs were written with dups of finalHdr).

  # Build repProv: compositeOff (representative) → ClusterProv.
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

  # chrom → header-order index for sort.
  var chromIdx: Table[string, int]
  for i, c in chromOrder.pairs: chromIdx[c] = i

  # Output-time INFO filter: drops always-keep-but-not-user-requested fields.
  # Provenance (SOURCE/SOURCELIST/N_SOURCE/N_MERGED) is always kept; MATCHA_BOFF
  # is always dropped.
  let infoFilter = cfg.infoFields
  proc keepInfoOut(name: string): bool =
    if name == "MATCHA_BOFF" or name == "MATCHA_CALLER_IDX": return false
    if name in ["SOURCE", "SOURCELIST", "N_SOURCE", "N_MERGED"]: return true
    if infoFilter.len == 0: return true
    for tok in infoFilter:
      if name == tok or name.startsWith(tok & "_"): return true
    false

  # Stream merged BCFs, collect representatives.
  var buf: seq[BufferedRep]
  var boffScratch: seq[int32]
  for path in mergedPaths.values:
    var vcf: VCF
    if not open(vcf, path):
      raise newException(IOError, "cannot open merged BCF: " & path)

    for v in vcf:
      let off = readBoff(v, boffScratch)
      if off notin repProv: continue
      let prov = repProv[off]

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
    if cfg.outputPath.endsWith(".bcf"):    "wb"
    elif cfg.outputPath.endsWith(".vcf.gz"): "wz"
    else:                                   "w"

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
      stderr.writeLine "error: input file not found: " & caller.path
      quit(1)
  if cfg.tmpDir == "":
    stderr.writeLine "error: tmpDir must be set"
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
        stderr.writeLine "error: cannot open: " & caller.path
        quit(1)
      counts.add(bcf_hdr_nsamples(vcf.header.hdr).int)
      vcf.close()
    for i, n in counts.pairs:
      if n > 1:
        stderr.writeLine "error: caller '" & cfg.callers[i].name & "' (" &
                         cfg.callers[i].path & ") has " & $n &
                         " sample columns; matcha collapse supports at most " &
                         "1 sample per input (split multi-sample VCFs first)"
        quit(1)
    let first = counts[0]
    for i, n in counts.pairs:
      if n != first:
        stderr.writeLine "error: inconsistent sample counts across inputs: " &
                         "caller '" & cfg.callers[0].name & "' has " & $first &
                         " sample(s) but caller '" & cfg.callers[i].name &
                         "' has " & $n &
                         "; all collapse inputs must have the same sample count"
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

  # Phase 2: integrated preproc+merge — stream all N caller VCFs in lockstep.
  logV("integrated preproc+merge over " & $cfg.callers.len & " caller(s)")
  let im = integratedMerge(cfg.callers, mh, chromOrder, cfg,
                           CollapseVersion, cmdLine)

  # Build a PreprocOutput describing the merged slim BCFs for buildWorkQueue.
  let mergedPreproc = PreprocOutput(
    paths:          im.paths,
    populatedBins:  im.populatedBins,
    chromsBySvtype: im.chromsBySvtype,
    chromOrder:     chromOrder,
  )

  # Phase 3: Pass 1 — self-match over merged slim BCFs.
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

  # Phase 4: Pass 2 — enumerate allOffsets + PASS/QUAL from merged BCFs.
  let simMap = buildSimilarityMap(allPairs)
  let (allOffsets, passQualMap) = exploreMerged(im.paths)
  logV("offsets seen: " & $allOffsets.len)

  # Phase 5: cluster and select representatives.
  let clusters = clusterAll(allOffsets, simMap, cfg.linkage, cfg.threshold)
  var finalClusters: seq[seq[int64]]
  for cl in clusters:
    let rep = selectRepresentative(cl, simMap, passQualMap, cfg.priority)
    var ordered = @[rep]
    for off in cl:
      if off != rep: ordered.add(off)
    finalClusters.add(ordered)
  logV("clusters: " & $finalClusters.len)

  # Phase 6: write output.
  writeOutput(cfg, im.finalHdr, chromOrder, im.paths, finalClusters)

  # Clean up: merged slim BCFs + CSI indexes, and finalHdr.
  for path in im.paths.values:
    if fileExists(path):          removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")
  bcf_hdr_destroy(im.finalHdr)
