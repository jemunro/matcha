## merge.nim — `matcha merge` subcommand.
##
## Cohort pVCF merge across N single-sample SV VCFs (typically `matcha
## collapse` outputs). Clusters equivalent SVs across samples, emits one
## row per site with N FORMAT columns, and computes cohort INFO (AC/AN/AF).
##
## Pipeline:
##   1. validate inputs (1 sample per file, distinct sample IDs)
##   2. resolveHeaders + buildSlimHdr + buildOutputHdr
##   3. integratedMerge (stamps FORMAT/SID, preserves BND ALT)
##   4. self-match (emitSingletons=true) over the merged slim BCFs
##   5. cluster + select representative (reused from collapse)
##   6. writeMergeOutput — fresh N-sample bcf1_t per cluster, AC/AN/AF
##      computed from assembled GT, CALLERS union when input had them.

import std/[algorithm, math, os, sequtils, sets, strutils, tables]
import hts
import hts/private/hts_concat
import utils, preproc, matchcore, mergecore, log, synced_bcf_reader

# ---------------------------------------------------------------------------
# MergeConfig
# ---------------------------------------------------------------------------

type
  MergeConfig* = object
    metric*:       Metric
    threshold*:    float64
    bndSlop*:      int
    insSlop*:      int
    insMinSim*:    float64
    linkage*:      LinkageMethod
    priority*:     seq[PriorityCriterion]
    formatFields*: seq[string]   ## FORMAT fields to carry; GT auto-added
    infoFields*:   seq[string]   ## --info filter; empty = auto-extracted only
    outputPath*:   string
    nThreads*:     int
    tmpDir*:       string
    callers*:      seq[CallerInput]
    keptChrs*:     seq[string]    ## --chrs filter; empty = no filter.
    missingToRef*: bool           ## --missing-to-ref: absent samples → 0/0.

# ---------------------------------------------------------------------------
# Sample-ID validation
# ---------------------------------------------------------------------------

proc validateMergeInputs(callers: seq[CallerInput]): seq[string] =
  ## Open each input, enforce exactly 1 sample, return sample IDs in CLI
  ## order. Rejects duplicate sample IDs.
  result = newSeq[string](callers.len)
  var seen: HashSet[string]
  for ci, caller in callers.pairs:
    var vcf: VCF
    if not open(vcf, caller.path):
      logError("cannot open: " & caller.path); quit(1)
    let n = bcf_hdr_nsamples(vcf.header.hdr).int
    if n != 1:
      logError("merge: input '" & caller.path & "' has " & $n &
               " sample column(s); merge requires exactly 1 per input")
      vcf.close(); quit(1)
    let s = $cast[cstringArray](vcf.header.hdr.samples)[0]
    if s in seen:
      logError("merge: duplicate sample ID '" & s & "' in input '" &
               caller.path & "'")
      vcf.close(); quit(1)
    seen.incl(s)
    result[ci] = s
    vcf.close()

# ---------------------------------------------------------------------------
# Header helpers
# ---------------------------------------------------------------------------

proc fmtFieldType(h: ptr bcf_hdr_t; name: string): cint =
  ## Returns BCF_HT_INT/REAL/STR for FORMAT field `name`, or -1 if absent.
  let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](h.hrec)
  for i in 0 ..< h.nhrec.int:
    let hr = hrecs[i]
    if hr.`type` != BCF_HEADER_TYPE.BCF_HL_FMT.cint: continue
    let keys = cast[ptr UncheckedArray[cstring]](hr.keys)
    let vals = cast[ptr UncheckedArray[cstring]](hr.vals)
    var id, typ: string
    for j in 0 ..< hr.nkeys.int:
      case $keys[j]
      of "ID":   id  = $vals[j]
      of "Type": typ = $vals[j]
    if id == name:
      case typ
      of "Integer", "Flag": return BCF_HT_INT.cint
      of "Float":           return BCF_HT_REAL.cint
      of "String", "Character": return BCF_HT_STR.cint
      else: return -1.cint
  -1.cint

# ---------------------------------------------------------------------------
# buildSlimHdr — header for per-(svtype,bin) merged slim BCFs
# ---------------------------------------------------------------------------

proc keepInfoForMergeOut(name: string; infoFilter: seq[string]): bool =
  ## Output INFO keep-set: drop SRC_INDEX/CALLER_IDX; always keep
  ## SVTYPE/SVLEN/END/CHR2/POS2 + cohort fields + CALLERS/N_CALLERS;
  ## filter the rest through --info.
  if name in ["SRC_INDEX", "CALLER_IDX"]: return false
  if name in ["SVTYPE", "SVLEN", "END", "CHR2", "POS2",
              "AC", "AN", "AF", "CALLERS", "N_CALLERS"]: return true
  if infoFilter.len == 0: return false
  for tok in infoFilter:
    if name == tok or name.startsWith(tok & "_"): return true
  false

proc buildSlimHdr(callers: seq[CallerInput]; mh: MergedHeader;
                  cfg: MergeConfig): ptr bcf_hdr_t =
  ## 1 dummy sample "SAMPLE" + user --format + FORMAT/SID +
  ## matchcore-required INFO + user --info.
  result = bcf_hdr_init("w".cstring)

  addContigsUnion(result, callers, toHashSet(cfg.keptChrs))
  addFiltersUnion(result, callers)

  # INFO from MergedHeader, filtered through keepInfoForMerged. FORMAT
  # restricted to user --format.
  let fmtKeep = toHashSet(cfg.formatFields)
  let userInfo = cfg.infoFields
  addHeaderLinesFiltered(result, mh,
    keepInfo = proc (id: string): bool =
      keepInfoForMerged(id, userInfo),
    keepFmt  = proc (id: string): bool =
      id in fmtKeep)

  addStandardSvInfoDefs(result)

  # Matcha-internal INFO defs.
  discard bcf_hdr_append(result,
    "##INFO=<ID=SRC_INDEX,Number=1,Type=Integer,Description=\"matcha-internal: sequential record index\">".cstring)
  discard bcf_hdr_append(result,
    "##INFO=<ID=CALLER_IDX,Number=1,Type=Integer,Description=\"matcha-internal: caller index (0-based)\">".cstring)

  # FORMAT/SID — slim-internal source sample ID.
  discard bcf_hdr_append(result,
    "##FORMAT=<ID=SID,Number=1,Type=String,Description=\"matcha-internal: source sample ID\">".cstring)
  # GT auto-added so the slim hdr can carry it even when the user didn't
  # request it explicitly (cohort INFO needs it).
  if "GT" notin fmtKeep:
    discard bcf_hdr_append(result,
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".cstring)

  # 1 dummy sample.
  discard bcf_hdr_add_sample(result, "SAMPLE".cstring)
  discard bcf_hdr_add_sample(result, cstring(nil))
  discard bcf_hdr_sync(result)

# ---------------------------------------------------------------------------
# buildOutputHdr — header for the final N-sample cohort pVCF
# ---------------------------------------------------------------------------

proc buildOutputHdr(callers: seq[CallerInput]; sampleIds: seq[string];
                    mh: MergedHeader; cfg: MergeConfig;
                    inputsHadCallers: bool;
                    cmdLine: string): ptr bcf_hdr_t =
  result = bcf_hdr_init("w".cstring)

  addContigsUnion(result, callers, toHashSet(cfg.keptChrs))
  addFiltersUnion(result, callers)

  # INFO from MergedHeader, filtered through keepInfoForMergeOut. FORMAT
  # restricted to user --format (no SID — slim-internal).
  let fmtKeep = toHashSet(cfg.formatFields)
  let userInfo = cfg.infoFields
  addHeaderLinesFiltered(result, mh,
    keepInfo = proc (id: string): bool =
      keepInfoForMergeOut(id, userInfo),
    keepFmt  = proc (id: string): bool =
      id in fmtKeep)

  addStandardSvInfoDefs(result)

  # Cohort INFO defs.
  let haveInfo = infoFieldIds(result)
  if "AC" notin haveInfo:
    discard bcf_hdr_append(result,
      "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alt allele count\">".cstring)
  if "AN" notin haveInfo:
    discard bcf_hdr_append(result,
      "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total alleles called\">".cstring)
  if "AF" notin haveInfo:
    discard bcf_hdr_append(result,
      "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alt allele frequency\">".cstring)

  # CALLERS / N_CALLERS (only when input records carried them).
  if inputsHadCallers:
    if "CALLERS" notin haveInfo:
      discard bcf_hdr_append(result,
        "##INFO=<ID=CALLERS,Number=.,Type=String,Description=\"Caller names contributing to this site (union across samples)\">".cstring)
    if "N_CALLERS" notin haveInfo:
      discard bcf_hdr_append(result,
        "##INFO=<ID=N_CALLERS,Number=1,Type=Integer,Description=\"Distinct callers contributing to this site\">".cstring)

  # GT def — auto-add if user --format omits it (cohort INFO needs it).
  if "GT" notin fmtFieldIds(result):
    discard bcf_hdr_append(result,
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".cstring)

  discard bcf_hdr_append(result,
    ("##source=matcha merge " & MatchaVersion).cstring)
  if cmdLine.len > 0:
    discard bcf_hdr_append(result, ("##matcha_cmdline=" & cmdLine).cstring)

  # N samples in CLI order.
  for s in sampleIds:
    discard bcf_hdr_add_sample(result, s.cstring)
  discard bcf_hdr_add_sample(result, cstring(nil))
  discard bcf_hdr_sync(result)

# ---------------------------------------------------------------------------
# Output assembly
# ---------------------------------------------------------------------------

type
  MemberRec = object
    rec:        ptr bcf1_t   ## owned dup of the slim record
    sid:        string
    callerIdx:  int32
    hasPASS:    bool
    qual:       uint16

  BufferedRow = object
    chromOrderIdx: int
    pos:           int64
    rec:           ptr bcf1_t

# Pick the highest-priority member when multiple records share a SID
# (rare same-sample cluster collision).
proc pickPriorityMember(members: seq[int]; recs: seq[MemberRec];
                        priority: seq[PriorityCriterion]): int =
  if members.len == 1: return members[0]
  var cands = members
  for crit in priority:
    if cands.len == 1: break
    case crit
    of pcPass:
      let p = cands.filterIt(recs[it].hasPASS)
      if p.len > 0: cands = p
    of pcQual:
      var bestQ = recs[cands[0]].qual
      for i in cands:
        if recs[i].qual > bestQ: bestQ = recs[i].qual
      cands = cands.filterIt(recs[it].qual >= bestQ)
    of pcCentre, pcOrder:
      var bestCi = high(int32)
      for i in cands:
        if recs[i].callerIdx < bestCi: bestCi = recs[i].callerIdx
      cands = cands.filterIt(recs[it].callerIdx == bestCi)
  cands[0]

# Reconstruct symbolic ALT for interval SVTYPEs.
proc symbolicAlt(svt: SvType): string =
  case svt
  of svDEL: "<DEL>"
  of svDUP: "<DUP>"
  of svINV: "<INV>"
  of svINS: "<INS>"
  else: "."

# Copy a single INFO field from src record (in srcHdr space) to dst record
# (in dstHdr space). Skips if absent on src.
proc copyInfoField(dstHdr, srcHdr: ptr bcf_hdr_t; dst, src: ptr bcf1_t;
                   name: string; buf: var pointer; bufN: var cint) =
  let info = bcf_get_info(srcHdr, src, name.cstring)
  if info == nil: return
  let htype =
    case info.`type`.cint
    of BCF_BT_FLOAT: BCF_HT_REAL.cint
    of BCF_BT_CHAR:  BCF_HT_STR.cint
    else:            BCF_HT_INT.cint
  let n = bcf_get_info_values(srcHdr, src, name.cstring,
                              buf.addr, bufN.addr, htype)
  if n > 0:
    discard bcf_update_info(dstHdr, dst, name.cstring, buf, n.cint, htype)

# ---------------------------------------------------------------------------
# writeMergeOutput
# ---------------------------------------------------------------------------

proc writeMergeOutput(cfg: MergeConfig;
                     outputHdr, slimHdrShared: ptr bcf_hdr_t;
                     sampleIds: seq[string];
                     sampleIdxBySID: Table[string, int];
                     chromOrder: seq[string];
                     mergedPaths: Table[SvtypeBin, string];
                     fileList: seq[string];
                     locByIdx: Table[int32, tuple[chromIdx: int16; pos: int32; fileIdx: int16]];
                     passQualMap: Table[int32, tuple[hasPASS: bool; qual: uint16; callerIdx: int32]];
                     finalClusters: seq[seq[int32]];
                     emitCallers: bool) =

  let nSamples = sampleIds.len
  let formatFields = cfg.formatFields

  # Pre-resolve INFO field types from outputHdr for copying.
  proc keepInfoOut(name: string): bool =
    keepInfoForMergeOut(name, cfg.infoFields)

  # Group cluster members by fileIdx so we open each slim BCF once.
  var membersByFile: Table[int16, seq[int32]]
  for cl in finalClusters:
    for idx in cl:
      let loc = locByIdx[idx]
      membersByFile.mgetOrPut(loc.fileIdx, @[]).add(idx)
  # Build a per-cluster member list of MemberRec, fetching records via CSI.
  var clusterMembers = newSeq[seq[int]](finalClusters.len)
  var allRecs: seq[MemberRec]
  # Build idx → clusterIdx lookup
  var idxToCluster: Table[int32, int]
  for ci, cl in finalClusters.pairs:
    for idx in cl:
      idxToCluster[idx] = ci

  # Per-fileIdx CSI sweep.
  var idxBuf, ciBuf: seq[int32]
  var sidBuf: pointer = nil
  var sidBufN: cint = 0
  for fileIdx, idxs in membersByFile.pairs:
    var vcf: VCF
    if not open(vcf, fileList[fileIdx]):
      raise newException(IOError, "cannot open merged BCF: " & fileList[fileIdx])
    # Sort idxs by pos to allow forward CSI streaming.
    var sorted = idxs
    sorted.sort(proc(a, b: int32): int =
      cmp(locByIdx[a].pos, locByIdx[b].pos))
    for idx in sorted:
      let loc = locByIdx[idx]
      let region = chromOrder[loc.chromIdx] & ":" & $loc.pos & "-" & $loc.pos
      var found = false
      for v in vcf.query(region):
        if readSrcIndex(v, idxBuf) != idx: continue
        # Read FORMAT/SID. `n` is the total byte count; for our 1-sample
        # slim records the entire buffer is that sample's string bytes
        # (possibly null-padded). Strip trailing nulls.
        let n = bcf_get_format_values(vcf.header.hdr, v.c, "SID".cstring,
                                      sidBuf.addr, sidBufN.addr,
                                      BCF_HT_STR.cint)
        var sid: string
        if n > 0 and sidBuf != nil:
          let raw = cast[ptr UncheckedArray[char]](sidBuf)
          var sLen = n.int
          while sLen > 0 and raw[sLen - 1] == '\0': dec sLen
          sid = newString(sLen)
          for k in 0 ..< sLen: sid[k] = raw[k]
        let ciVal =
          if v.info().get("CALLER_IDX", ciBuf) == Status.OK and ciBuf.len > 0:
            ciBuf[0]
          else: 0'i32
        let pq = passQualMap.getOrDefault(idx,
                   (hasPASS: false, qual: 0'u16, callerIdx: ciVal))
        let dup = bcf_dup(v.c)
        let memIdx = allRecs.len
        allRecs.add(MemberRec(rec: dup, sid: sid, callerIdx: pq.callerIdx,
                              hasPASS: pq.hasPASS, qual: pq.qual))
        clusterMembers[idxToCluster[idx]].add(memIdx)
        found = true
        break
      if not found:
        logWarn("merge: could not retrieve slim record for SRC_INDEX " & $idx)
    vcf.close()
  if sidBuf != nil: c_free(sidBuf)

  # chrom name → output-order index for sorting.
  var chromIdx: Table[string, int]
  for i, c in chromOrder.pairs: chromIdx[c] = i

  # Buffer output records, then sort+write.
  var bufRows: seq[BufferedRow]

  # Scratch buffers reused across clusters.
  var infoCopyBuf: pointer = nil
  var infoCopyN: cint = 0
  var fmtReadBuf: pointer = nil
  var fmtReadN: cint = 0

  for clusterIdx, cl in finalClusters.pairs:
    if cl.len == 0: continue
    let repIdx = cl[0]
    let members = clusterMembers[clusterIdx]
    if members.len == 0:
      logWarn("merge: cluster " & $clusterIdx & " has no retrieved members; skipping")
      continue

    # Find the MemberRec for the representative.
    var repMemberPos = -1
    for k, mi in members.pairs:
      # repIdx is the SRC_INDEX; we need to identify which member corresponds.
      # idxToCluster keys are SRC_INDEX → cluster, but we don't have the
      # reverse here without a lookup. Re-establish: scan members for the
      # one whose dup record has SRC_INDEX == repIdx.
      var idxBufLocal: seq[int32]
      var v: Variant
      # We can extract SRC_INDEX from the dup using a Variant wrapper —
      # cheaper: store srcIndex on MemberRec. (Patched below.)
      discard v
      discard idxBufLocal
      if allRecs[mi].callerIdx >= 0:
        # Defer the actual rep-identification to a second pass; for now
        # mark and patch below.
        discard
    # Linear scan with a SRC_INDEX extraction (no scratch on MemberRec).
    block findRep:
      var sb: seq[int32]
      for mi in members:
        # Read SRC_INDEX off the duped record.
        var nval: cint = 0
        var p: pointer = nil
        let nr = bcf_get_info_values(slimHdrShared, allRecs[mi].rec,
                                     "SRC_INDEX".cstring,
                                     p.addr, nval.addr, BCF_HT_INT.cint)
        if nr >= 1 and p != nil:
          let arr = cast[ptr UncheckedArray[int32]](p)
          let s = arr[0]
          if s == repIdx:
            repMemberPos = mi
            c_free(p)
            break findRep
          c_free(p)
        elif p != nil:
          c_free(p)
        discard sb
    if repMemberPos < 0:
      # Fallback: first member.
      repMemberPos = members[0]

    let repRec = allRecs[repMemberPos].rec

    # --- Assemble slot map: SID → memberRec index.
    var slotByMember = newSeq[int](nSamples)
    for i in 0 ..< nSamples: slotByMember[i] = -1
    var bySid: Table[string, seq[int]]
    for mi in members:
      if allRecs[mi].sid.len == 0: continue
      bySid.mgetOrPut(allRecs[mi].sid, @[]).add(mi)
    for sid, list in bySid.pairs:
      let s = sampleIdxBySID.getOrDefault(sid, -1)
      if s < 0:
        logWarn("merge: cluster " & $clusterIdx &
                ": record SID '" & sid & "' has no output sample column")
        continue
      let pick =
        if list.len == 1: list[0]
        else:
          pickPriorityMember(list, allRecs, cfg.priority)
      slotByMember[s] = pick

    # --- Build fresh output record.
    let outRec = bcf_init()
    outRec.rid = repRec.rid
    outRec.pos = repRec.pos
    outRec.qual = repRec.qual

    # ID — copy from rep slim record.
    if repRec.d.id != nil:
      discard bcf_update_id(outputHdr, outRec, repRec.d.id)

    # SVTYPE — read from rep INFO (needed for symbolic ALT reconstruction).
    var svtBuf: pointer = nil
    var svtN: cint = 0
    let svtRead = bcf_get_info_values(slimHdrShared, repRec, "SVTYPE".cstring,
                                      svtBuf.addr, svtN.addr, BCF_HT_STR.cint)
    var svt: SvType = svUNKNOWN
    if svtRead > 0 and svtBuf != nil:
      svt = parseSvType($cast[cstring](svtBuf))
    if svtBuf != nil: c_free(svtBuf)

    # REF/ALT. For BND and sequence-resolved INS the slim record preserved
    # the source REF/ALT verbatim; emit those. For other types (and INS that
    # was symbolic on input) fall back to a symbolic ALT.
    if svt == svBND or svt == svINS:
      discard bcf_unpack(repRec, BCF_UN_STR.cint)
      let alleles = cast[ptr UncheckedArray[cstring]](repRec.d.allele)
      let preserved = repRec.n_allele >= 2 and
                      (alleles == nil or alleles[1] == nil or
                       ($alleles[1]).len > 0)
      if preserved and repRec.n_allele >= 2:
        var altStr = $alleles[0]
        for i in 1 ..< repRec.n_allele.int:
          altStr &= ","
          altStr &= $alleles[i]
        discard bcf_update_alleles_str(outputHdr, outRec, altStr.cstring)
      elif svt == svINS:
        discard bcf_update_alleles_str(outputHdr, outRec, "N,<INS>".cstring)
      else:
        discard bcf_update_alleles_str(outputHdr, outRec, "N,.".cstring)
    else:
      discard bcf_update_alleles_str(outputHdr, outRec,
        ("N," & symbolicAlt(svt)).cstring)

    # FILTER — translate IDs from slim hdr space to outputHdr space.
    discard bcf_unpack(repRec, BCF_UN_FLT.cint)
    if repRec.d.n_flt > 0:
      let fltArr = cast[ptr UncheckedArray[cint]](repRec.d.flt)
      for j in 0 ..< repRec.d.n_flt.int:
        let name = $hdrInt2Id(slimHdrShared, BCF_DT_ID.cint, fltArr[j])
        let newId = bcf_hdr_id2int(outputHdr, BCF_DT_ID.cint, name.cstring)
        if newId >= 0:
          discard bcf_add_filter(outputHdr, outRec, newId)

    # INFO — copy non-internal fields from rep through keepInfoOut filter.
    discard bcf_unpack(repRec, BCF_UN_INFO.cint)
    let idPairs = cast[ptr UncheckedArray[bcf_idpair_t]](slimHdrShared.id[BCF_DT_ID])
    if repRec.n_info > 0:
      let infoArr = cast[ptr UncheckedArray[bcf_info_t]](repRec.d.info)
      for i in 0 ..< repRec.n_info.int:
        let key = infoArr[i].key
        if key < 0: continue
        let nm = $idPairs[key].key
        if not keepInfoOut(nm): continue
        if nm in ["AC", "AN", "AF"]: continue  # computed below
        if nm in ["CALLERS", "N_CALLERS"]: continue  # unioned below
        copyInfoField(outputHdr, slimHdrShared, outRec, repRec, nm,
                      infoCopyBuf, infoCopyN)

    # --- Build per-sample FORMAT data field-by-field, computing AC/AN from GT.
    var anTotal = 0
    var acTotal = 0
    for fieldName in formatFields:
      # GT is always encoded as Int regardless of header Type=String.
      let htype =
        if fieldName == "GT": BCF_HT_INT.cint
        else: fmtFieldType(outputHdr, fieldName)
      if htype < 0:
        # Field not declared in output header — skip.
        continue

      # Read each slot's FORMAT data; determine maxK; assemble N*maxK.
      var perSampleInt: seq[seq[int32]] = newSeq[seq[int32]](nSamples)
      var perSampleFloat: seq[seq[float32]] = newSeq[seq[float32]](nSamples)
      var perSampleStr: seq[string] = newSeq[string](nSamples)
      var maxK = 0
      for s in 0 ..< nSamples:
        let mi = slotByMember[s]
        if mi < 0: continue
        let memRec = allRecs[mi].rec
        case htype.int
        of BCF_HT_INT.int:
          fmtReadN = 0; fmtReadBuf = nil
          let n = bcf_get_format_values(slimHdrShared, memRec,
                                        fieldName.cstring,
                                        fmtReadBuf.addr, fmtReadN.addr,
                                        BCF_HT_INT.cint)
          if n > 0 and fmtReadBuf != nil:
            let arr = cast[ptr UncheckedArray[int32]](fmtReadBuf)
            var vals = newSeq[int32](n.int)
            for k in 0 ..< n.int: vals[k] = arr[k]
            perSampleInt[s] = vals
            if vals.len > maxK: maxK = vals.len
          if fmtReadBuf != nil:
            c_free(fmtReadBuf); fmtReadBuf = nil; fmtReadN = 0
        of BCF_HT_REAL.int:
          fmtReadN = 0; fmtReadBuf = nil
          let n = bcf_get_format_values(slimHdrShared, memRec,
                                        fieldName.cstring,
                                        fmtReadBuf.addr, fmtReadN.addr,
                                        BCF_HT_REAL.cint)
          if n > 0 and fmtReadBuf != nil:
            let arr = cast[ptr UncheckedArray[float32]](fmtReadBuf)
            var vals = newSeq[float32](n.int)
            for k in 0 ..< n.int: vals[k] = arr[k]
            perSampleFloat[s] = vals
            if vals.len > maxK: maxK = vals.len
          if fmtReadBuf != nil:
            c_free(fmtReadBuf); fmtReadBuf = nil; fmtReadN = 0
        of BCF_HT_STR.int:
          fmtReadN = 0; fmtReadBuf = nil
          let n = bcf_get_format_values(slimHdrShared, memRec,
                                        fieldName.cstring,
                                        fmtReadBuf.addr, fmtReadN.addr,
                                        BCF_HT_STR.cint)
          if n > 0 and fmtReadBuf != nil:
            perSampleStr[s] = $cast[cstring](fmtReadBuf)
          if fmtReadBuf != nil:
            c_free(fmtReadBuf); fmtReadBuf = nil; fmtReadN = 0
        else: discard

      # Write N-sample buffer.
      case htype.int
      of BCF_HT_INT.int:
        if maxK == 0: maxK = 1
        # GT uses bcf_gt_missing (= 0) for missing alleles, not INT32_MIN.
        # With --missing-to-ref, absent samples get encoded REF (= 2,
        # i.e. (0+1)<<1 | 0 = REF allele unphased) instead.
        let missingVal: int32 =
          if fieldName == "GT":
            if cfg.missingToRef: 2'i32 else: 0'i32
          else: bcfInt32Missing
        var buf = newSeq[int32](nSamples * maxK)
        for s in 0 ..< nSamples:
          let mi = slotByMember[s]
          if mi < 0:
            # Missing sample.
            for k in 0 ..< maxK:
              buf[s * maxK + k] = missingVal
          else:
            let vals = perSampleInt[s]
            for k in 0 ..< maxK:
              if k < vals.len:
                buf[s * maxK + k] = vals[k]
              else:
                buf[s * maxK + k] = bcfInt32VectorEnd
        discard bcf_update_format(outputHdr, outRec, fieldName.cstring,
                                  buf[0].addr, (nSamples * maxK).cint,
                                  BCF_HT_INT.cint)
        # GT contribution to AC/AN.
        if fieldName == "GT":
          for s in 0 ..< nSamples:
            for k in 0 ..< maxK:
              let v = buf[s * maxK + k]
              if v == bcfInt32VectorEnd: continue
              if v == 0: continue
              let allele = (v shr 1) - 1
              inc anTotal
              if allele >= 1: inc acTotal
      of BCF_HT_REAL.int:
        if maxK == 0: maxK = 1
        var buf = newSeq[float32](nSamples * maxK)
        for s in 0 ..< nSamples:
          let mi = slotByMember[s]
          if mi < 0:
            for k in 0 ..< maxK:
              buf[s * maxK + k] = bcfFloatMissing()
          else:
            let vals = perSampleFloat[s]
            for k in 0 ..< maxK:
              if k < vals.len:
                buf[s * maxK + k] = vals[k]
              else:
                buf[s * maxK + k] = bcfFloatVectorEnd()
        discard bcf_update_format(outputHdr, outRec, fieldName.cstring,
                                  buf[0].addr, (nSamples * maxK).cint,
                                  BCF_HT_REAL.cint)
      of BCF_HT_STR.int:
        # cstringArray of N values (missing → ".").
        var strs = newSeq[string](nSamples)
        var cstrs = newSeq[cstring](nSamples)
        for s in 0 ..< nSamples:
          let mi = slotByMember[s]
          strs[s] = if mi < 0: "." else: perSampleStr[s]
          cstrs[s] = strs[s].cstring
        discard bcf_update_format_string(outputHdr, outRec,
                                         fieldName.cstring,
                                         cast[cstringArray](cstrs[0].addr),
                                         nSamples.cint)
      else: discard

    # --- Cohort INFO: AC/AN/AF.
    var anVal = anTotal.int32
    var acVal = acTotal.int32
    discard bcf_update_info(outputHdr, outRec, "AN".cstring,
                            anVal.addr, 1.cint, BCF_HT_INT.cint)
    discard bcf_update_info(outputHdr, outRec, "AC".cstring,
                            acVal.addr, 1.cint, BCF_HT_INT.cint)
    if anTotal > 0:
      var afVal = float32(acTotal) / float32(anTotal)
      discard bcf_update_info(outputHdr, outRec, "AF".cstring,
                              afVal.addr, 1.cint, BCF_HT_REAL.cint)
    else:
      var afMissing = bcfFloatMissing()
      discard bcf_update_info(outputHdr, outRec, "AF".cstring,
                              afMissing.addr, 1.cint, BCF_HT_REAL.cint)

    # --- CALLERS union (only when inputs carried CALLERS in their INFO).
    if emitCallers:
      var unionList: seq[string]
      var seen: HashSet[string]
      # Representative first, then rest in caller-idx order.
      var ordered: seq[int] = @[repMemberPos]
      for mi in members:
        if mi != repMemberPos: ordered.add(mi)
      var callerBuf: pointer = nil
      var callerBufN: cint = 0
      for mi in ordered:
        let n = bcf_get_info_values(slimHdrShared, allRecs[mi].rec,
                                    "CALLERS".cstring,
                                    callerBuf.addr, callerBufN.addr,
                                    BCF_HT_STR.cint)
        if n > 0 and callerBuf != nil:
          let raw = $cast[cstring](callerBuf)
          for tok in raw.split(','):
            let t = tok.strip
            if t.len > 0 and t notin seen:
              seen.incl(t)
              unionList.add(t)
        if callerBuf != nil:
          c_free(callerBuf); callerBuf = nil; callerBufN = 0
      if unionList.len > 0:
        var joined = unionList.join(",")
        discard bcf_update_info(outputHdr, outRec, "CALLERS".cstring,
                                joined[0].addr, joined.len.cint,
                                BCF_HT_STR.cint)
        var nc = unionList.len.int32
        discard bcf_update_info(outputHdr, outRec, "N_CALLERS".cstring,
                                nc.addr, 1.cint, BCF_HT_INT.cint)

    let chromName = $hdrInt2Id(outputHdr, BCF_DT_CTG.cint, outRec.rid)
    let coi = chromIdx.getOrDefault(chromName, high(int))
    bufRows.add(BufferedRow(chromOrderIdx: coi, pos: outRec.pos, rec: outRec))

  if infoCopyBuf != nil: c_free(infoCopyBuf)

  # Free per-cluster member dups.
  for mr in allRecs:
    if mr.rec != nil: bcf_destroy(mr.rec)

  # Sort + write.
  bufRows.sort(proc(a, b: BufferedRow): int =
    let c = cmp(a.chromOrderIdx, b.chromOrderIdx)
    if c != 0: c else: cmp(a.pos, b.pos))

  let outPath = if isStdoutPath(cfg.outputPath): "/dev/stdout" else: cfg.outputPath
  let mode =
    if cfg.outputPath.endsWith(".bcf"):       "wb"
    elif cfg.outputPath.endsWith(".vcf.gz"):  "wz"
    else:                                     "w"

  var outVcf: VCF
  if not open(outVcf, outPath, mode = mode):
    raise newException(IOError, "cannot open output: " & outPath)
  block:
    var dummy: VCF
    discard open(dummy, cfg.callers[0].path)
    outVcf.copy_header(dummy.header)
    dummy.close()
  bcf_hdr_destroy(outVcf.header.hdr)
  outVcf.header.hdr = outputHdr  # SHARED — clear before close to avoid double-free
  discard outVcf.write_header()

  # bgzf_mt is intentionally skipped when building an inline CSI index: in MT
  # mode bgzf_tell returns a stale block_address (not updated until the worker
  # thread flushes the block), corrupting all virtual offsets beyond the first
  # BGZF block.  Revisit using bcf_idx_init/bcf_idx_save for MT-safe indexing.
  let isBgzf = not isStdoutPath(cfg.outputPath) and
               (cfg.outputPath.endsWith(".bcf") or cfg.outputPath.endsWith(".vcf.gz"))
  var outIdx: ptr hts_idx_t = nil
  if isBgzf:
    let headerOff = uint64(bgzf_tell(bgzfHandle(outVcf)))
    outIdx = hts_idx_init(0, HTS_FMT_CSI.cint, headerOff, 14, 5)
    if outIdx == nil:
      raise newException(IOError, "cannot create CSI index for: " & outPath)

  for br in bufRows:
    discard bcf_write(vcfHtsFile(outVcf), outputHdr, br.rec)
    if outIdx != nil:
      let woff = uint64(bgzf_tell(bgzfHandle(outVcf)))
      discard hts_idx_push(outIdx, br.rec.rid, int64(br.rec.pos),
                           int64(br.rec.pos) + int64(br.rec.rlen), woff, 1)
    bcf_destroy(br.rec)

  if outIdx != nil:
    let finalOff = uint64(bgzf_tell(bgzfHandle(outVcf)))
    hts_idx_finish(outIdx, finalOff)
  outVcf.header.hdr = nil
  outVcf.close()
  logInfo("merge: wrote " & $bufRows.len & " record(s) to " &
       (if isStdoutPath(cfg.outputPath): "stdout" else: cfg.outputPath))
  if outIdx != nil:
    hts_idx_save(outIdx, cfg.outputPath.cstring, HTS_FMT_CSI.cint)
    hts_idx_destroy(outIdx)
    logInfo("indexed " & cfg.outputPath)

# ---------------------------------------------------------------------------
# runMerge — top-level entry point
# ---------------------------------------------------------------------------

proc runMerge*(cfg: MergeConfig; cmdLine: string = "") =
  if cfg.callers.len < 2:
    logError("merge: need at least 2 input files (got " & $cfg.callers.len & ")")
    quit(1)

  logInfo("matcha merge: " & $cfg.callers.len & " sample(s)" &
       " linkage=" & $cfg.linkage & " threads=" & $cfg.nThreads)

  for caller in cfg.callers:
    if not fileExists(caller.path):
      logError("input file not found: " & caller.path); quit(1)
    if not fileExists(caller.path & ".csi") and
       not fileExists(caller.path & ".tbi"):
      logError("no index found for: " & caller.path &
               " (run: bcftools index " & caller.path & ")"); quit(1)
  if cfg.tmpDir == "":
    logError("tmpDir must be set"); quit(1)

  # Sample validation + GT auto-add.
  let sampleIds = validateMergeInputs(cfg.callers)
  var cfgMut = cfg
  if "GT" notin cfgMut.formatFields:
    logInfo("merge: auto-adding GT to --format for cohort AC/AN/AF")
    cfgMut.formatFields = @["GT"] & cfgMut.formatFields

  # Resolve headers.
  let mh = resolveHeaders(cfgMut.callers,
    infoFilter = cfgMut.infoFields,
    fmtFilter  = cfgMut.formatFields)
  for w in mh.warnings: logWarn("merge header: " & w)

  # Detect whether any input declares INFO/CALLERS — drives output INFO defs.
  var inputsHadCallers = false
  for caller in cfgMut.callers:
    var vcf: VCF
    if not open(vcf, caller.path): continue
    if "CALLERS" in infoFieldIds(vcf.header.hdr): inputsHadCallers = true
    vcf.close()
    if inputsHadCallers: break

  # Chrom order anchored to first caller.
  var orderVcf: VCF
  if not open(orderVcf, cfgMut.callers[0].path):
    raise newException(IOError, "cannot open: " & cfgMut.callers[0].path)
  let chromOrder = captureChromOrder(orderVcf.header)
  orderVcf.close()

  # Build slim + output headers.
  let slimHdr   = buildSlimHdr(cfgMut.callers, mh, cfgMut)
  let outputHdr = buildOutputHdr(cfgMut.callers, sampleIds, mh, cfgMut,
                                  inputsHadCallers, cmdLine)

  # Integrated preproc + merge using slimHdr as writer.
  logInfo("integrated preproc+merge over " & $cfgMut.callers.len & " sample(s)")
  let msc = MergeStreamConfig(formatFields:     cfgMut.formatFields,
                               nThreads:         cfgMut.nThreads,
                               tmpDir:           cfgMut.tmpDir,
                               stampSID:         true,
                               sampleIdByCaller: sampleIds,
                               preserveBndAlt:   true,
                               preserveInsAlt:   true,
                               keptChrs:         toHashSet(cfgMut.keptChrs))
  let im = integratedMerge(cfgMut.callers, mh, slimHdr, msc, chromOrder)
  if cfgMut.keptChrs.len > 0:
    var seen: HashSet[string]
    for caller in cfgMut.callers:
      var v: VCF
      if not open(v, caller.path): continue
      for c in captureChromOrder(v.header): seen.incl(c)
      v.close()
    warnMissingChrs(cfgMut.keptChrs, seen)

  # Build a PreprocOutput describing the merged slim BCFs.
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

  # Self-match + cluster + representative selection (shared pipeline).
  logInfo("self-matching merged slim BCFs")
  let matchCfg = MatchConfig(
    metric:         cfgMut.metric,
    threshold:      cfgMut.threshold,
    bndSlop:        cfgMut.bndSlop,
    insSlop:        cfgMut.insSlop,
    insMinSim:      cfgMut.insMinSim,
    nThreads:       cfgMut.nThreads,
    tmpDir:         cfgMut.tmpDir,
    selfMode:       true,
    emitSingletons: true,
  )
  let cpr = selfMatchAndCluster(mergedPreproc, matchCfg,
                                 cfgMut.linkage, cfgMut.threshold,
                                 cfgMut.priority, "merge self-match",
                                 warnCallerStats = false)

  # SID → output sample index table.
  var sampleIdxBySID: Table[string, int]
  for i, s in sampleIds.pairs: sampleIdxBySID[s] = i

  # Output.
  writeMergeOutput(cfgMut, outputHdr, slimHdr, sampleIds, sampleIdxBySID,
                  chromOrder, im.paths, cpr.fileList, cpr.locByIdx,
                  cpr.passQualMap, cpr.finalClusters,
                  inputsHadCallers)

  # Teardown.
  removeDir(cfgMut.tmpDir)
  bcf_hdr_destroy(slimHdr)
  bcf_hdr_destroy(outputHdr)
