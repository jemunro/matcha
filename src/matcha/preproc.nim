## preproc.nim — VCF/BCF normalization, slimming, and per-(SVTYPE, bin) BCF output.
##
## For each input record we:
##   * Resolve SVTYPE (INFO/SVTYPE; symbolic ALT fallback; ALT wins on conflict).
##   * Resolve END   (INFO/END; POS+abs(SVLEN) fallback; require END > POS).
##   * Resolve SVLEN (INFO/SVLEN abs; END-POS fallback; warn+normalize on >10%
##     disagreement when both INFO fields were independently provided).
##   * Synthesize ID if absent ("." or empty) → CHROM_POS_SVTYPE_LINENUMBER.
##   * Slim INFO to the per-SVTYPE keep-set:
##     intervals → {END, SRC_INDEX}; BND → {CHR2, POS2, SRC_INDEX}.
##   * Compute size-bin via bins.binIndexFor(SVLEN).
##   * Write to a per-(SVTYPE, bin) BCF spanning all chroms; CSI-indexed.
##
## Skip categories (counted; final summary always emitted):
##   skUnsupportedSvtype  - INS only (silent count; out of scope)
##   skUnsupportedTra     - TRA (warn-emitting; not supported)
##   skUnresolvableSvtype - neither INFO/SVTYPE nor recognizable ALT
##   skUnresolvableEnd    - missing both INFO/END and INFO/SVLEN
##   skEndLePos           - resolved END <= POS
##   skMalformedBnd       - BND ALT could not be parsed for mate position
##   skUnknownContig      - BCF_ERR_CTG_UNDEF on read
##
## Per-record warnings are throttled at MATCHA_WARN_CAP (default 5) per reason.

import std/[algorithm, os, sequtils, sets, strutils, tables, tempfiles]
import hts
import hts/private/hts_concat  # BCF_ERR_CTG_UNDEF, BGZF, htsFile, bgzf_tell
import utils, log, bins, synced_bcf_reader

# hts-nim's VCF type keeps the underlying htsFile* private. We need access to
# it (specifically vcf.hts.fp.bgzf) so we can call bgzf_tell for CSI index
# building. Rather than patching vendored hts-nim, mirror enough of VCF's
# prefix to cast across. Both inherit RootObj, so the type-info slot lines up;
# the very next slot is `hts: ptr htsFile`.
# Layout is documented at vendor/hts-nim/src/hts/vcf.nim:19-26.
type
  VcfPriv = ref object of RootObj
    hts: ptr htsFile

proc bgzfHandle*(v: VCF): ptr BGZF {.inline.} =
  cast[VcfPriv](v).hts.fp.bgzf

proc vcfHtsFile*(v: VCF): ptr htsFile {.inline.} =
  cast[VcfPriv](v).hts

type
  SvtypeBin* = tuple[svtype: SvType, bin: int]

  PreprocOutput* = object
    paths*:          Table[SvtypeBin, string]         ## (svtype, bin) → BCF path
    populatedBins*:  Table[SvType, set[uint8]]        ## svtype → set of bin indexes
    chromsBySvtype*: Table[SvType, HashSet[string]]   ## svtype → chroms seen
    chromOrder*:     seq[string]                      ## chroms in input header order

  BinBEntry* = tuple[path: string; fileIdx: int16]

  MatchJob* = object
    chrom*:    string
    chromIdx*: int16                                   ## Index into chromOrder (for MatchPair.chromIdx).
    svtype*:   SvType
    binA*:     int
    pathA*:    string
    fileIdxA*: int16                                   ## Index into the run's slim-BCF file list.
    binsB*:    Table[int, BinBEntry]                   ## binB → (B's path, B's file index)

  SkipReason* = enum
    skUnsupportedSvtype
    skUnsupportedTra
    skUnresolvableSvtype
    skUnresolvableEnd
    skEndLePos
    skMalformedBnd
    skUnknownContig

  WarnState* = object
    skipped*:      array[SkipReason, int]
    emitted*:      array[SkipReason, int]
    syntheticId*:  int
    nRead*:        int
    nKept*:        int
    cap*:          int
    callset*:      string

proc initWarnState*(callset: string): WarnState =
  WarnState(cap: warnCap(), callset: callset)

# INFO fields kept in each class of temp BCF. SVTYPE is encoded in the filename;
# SVLEN is derivable from END-POS and matchcore never reads it.
const IntervalInfoDefs = [
  ("END",       "1", "Integer", "End position of the SV (1-based, inclusive)"),
  ("SRC_INDEX", "1", "Integer", "matcha-internal: sequential record index from preproc"),
]
const BndInfoDefs = [
  ("CHR2",      "1", "String",  "Chromosome of mate breakend"),
  ("POS2",      "1", "Integer", "Position of mate breakend"),
  ("SRC_INDEX", "1", "Integer", "matcha-internal: sequential record index from preproc"),
]

# Standard SV INFO definitions to backfill when the input header omits them.
# Without these, info.set() at write time fails because hts-nim resolves
# the field's type via the (reader's) header. Applied to the INPUT VCF only.
const SvInfoDefs = [
  ("SVTYPE",    "1", "String",  "Type of structural variant"),
  ("SVLEN",     "1", "Integer", "Length of the SV (absolute value, positive)"),
  ("END",       "1", "Integer", "End position of the SV (1-based, inclusive)"),
  ("CHR2",      "1", "String",  "Chromosome of mate breakend (BND/TRA)"),
  ("POS2",      "1", "Integer", "Position of mate breakend (BND/TRA)"),
  ("SRC_INDEX", "1", "Integer", "matcha-internal: sequential record index"),
]


proc reasonStr(r: SkipReason): string =
  case r
  of skUnsupportedSvtype:   "unsupported_svtype"
  of skUnsupportedTra:      "unsupported_tra"
  of skUnresolvableSvtype:  "unresolvable_svtype"
  of skUnresolvableEnd:     "unresolvable_end"
  of skEndLePos:            "end_le_pos"
  of skMalformedBnd:        "malformed_bnd"
  of skUnknownContig:       "unknown_contig"

proc warnSkip*(reason: SkipReason, ws: var WarnState,
              chrom: string, pos: int64, id: string, detail: string) =
  inc ws.skipped[reason]
  if ws.emitted[reason] < ws.cap:
    inc ws.emitted[reason]
    logWarn(chrom & ":" & $pos & " ID=" & id &
            " reason=\"" & reasonStr(reason) & ": " & detail & "\"")

proc symbolicAltSvtype(alt: string): SvType =
  ## Parse `<DEL>`, `<DUP>`, `<INV>`, `<INS>`, `<TRA>` → corresponding SvType.
  ## Returns svUNKNOWN for non-symbolic ALTs (e.g. BND breakend strings).
  if alt.len < 3 or alt[0] != '<' or alt[^1] != '>':
    return svUNKNOWN
  parseSvType(alt[1 .. ^2])

proc isBndAlt*(alt: string): bool =
  ## Detect VCF breakend ALT notation (`t[p[`, `t]p]`, `[p[t`, `]p]t`).
  ## Distinct from symbolic `<BND>` (handled by symbolicAltSvtype).
  if alt.len < 4 or alt[0] == '<': return false
  for c in alt:
    if c == '[' or c == ']': return true
  false

proc parseBndAlt*(alt: string): tuple[ok: bool; chr2: string; pos2: int64] =
  ## Extract (chr2, pos2) from a breakend ALT. Returns ok=false on any
  ## malformed input (missing bracket pair, missing colon, non-numeric pos,
  ## non-positive pos). Strand is intentionally ignored.
  if alt.len < 4: return
  var bracket: char
  var bStart = -1
  for i, c in alt:
    if c == '[' or c == ']':
      bracket = c; bStart = i; break
  if bStart < 0: return
  var bEnd = -1
  for i in bStart + 1 ..< alt.len:
    if alt[i] == bracket:
      bEnd = i; break
  if bEnd < 0 or bEnd <= bStart + 1: return
  let inner = alt[bStart + 1 ..< bEnd]
  let colon = inner.rfind(':')
  if colon < 1 or colon >= inner.len - 1: return
  let chrom = inner[0 ..< colon]
  let posStr = inner[colon + 1 .. ^1]
  if chrom.len == 0: return
  var pos: int64
  try:
    pos = parseBiggestInt(posStr).int64
  except ValueError:
    return
  if pos <= 0: return
  (ok: true, chr2: chrom, pos2: pos)

proc synthesizeId(chrom: string, pos: int64, svt: SvType, lineno: int): string =
  chrom & "_" & $pos & "_" & $svt & "_" & $lineno

# ---------------------------------------------------------------------------
# Raw C-API normalization helpers (used by preprocessVcf and integratedMerge)
# ---------------------------------------------------------------------------

proc getChromName*(hdr: ptr bcf_hdr_t; rid: int32): string =
  ## Look up contig name for `rid` directly from the header dictionary.
  if rid < 0: return ""
  let pairs = cast[ptr UncheckedArray[bcf_idpair_t]](hdr.id[BCF_DT_CTG])
  if pairs == nil: return ""
  result = $pairs[rid].key

proc getRecId*(rec: ptr bcf1_t): string =
  ## Return the record's ID (empty string for "." or missing).
  ## Caller must have called bcf_unpack(rec, BCF_UN_STR) first.
  if rec.d.id == nil: return ""
  let s = $rec.d.id
  if s == ".": "" else: s

proc resolveSvtype(v: Variant, infoBuf: var string): SvType =
  ## Resolve SVTYPE. ALT wins on conflict. Recognized ALT forms:
  ##   - symbolic `<DEL>` / `<DUP>` / `<INV>` / `<INS>` / `<TRA>` / `<BND>`
  ##   - bracket-form breakend (any non-symbolic ALT containing `[` or `]`)
  ## INFO/SVTYPE is the fallback. Returns svUNKNOWN if neither resolves.
  let infoOk = v.info().get("SVTYPE", infoBuf) == Status.OK and infoBuf.len > 0
  let infoSv =
    if infoOk: parseSvType($cast[cstring](addr infoBuf[0]))
    else:      svUNKNOWN
  var altSv = svUNKNOWN
  let alts = v.ALT
  if alts.len > 0:
    altSv = symbolicAltSvtype(alts[0])
    if altSv == svUNKNOWN and isBndAlt(alts[0]):
      altSv = svBND
  if altSv != svUNKNOWN:
    return altSv
  return infoSv

proc resolveEnd(v: Variant, infoEnd, infoSvlen: var seq[int32]):
    tuple[ok: bool, endPos: int64, fromInfo: bool] =
  if v.info().get("END", infoEnd) == Status.OK and infoEnd.len > 0:
    return (true, int64(infoEnd[0]), true)
  if v.info().get("SVLEN", infoSvlen) == Status.OK and infoSvlen.len > 0:
    return (true, v.POS + abs(int64(infoSvlen[0])), false)
  return (false, 0'i64, false)

proc resolveSvlen(v: Variant, endPos: int64, infoSvlen: var seq[int32]):
    tuple[ok: bool, svlen: int64, fromInfo: bool] =
  if v.info().get("SVLEN", infoSvlen) == Status.OK and infoSvlen.len > 0:
    return (true, abs(int64(infoSvlen[0])), true)
  let derived = endPos - v.POS
  if derived > 0:
    return (true, derived, false)
  return (false, 0'i64, false)

proc resolveRecord(view: Variant; ws: var WarnState;
                   chrom: string; pos1: int64; curId: string;
                   svtypeBuf: var string;
                   endBuf, svlenBuf: var seq[int32]):
    tuple[ok: bool; svt: SvType; endPos: int64; svlen: int64; binIdx: int;
          bndChr2: string; bndPos2: int64] =
  ## Resolve SVTYPE/END/SVLEN and skip-classify a record. Emits warnSkip for
  ## skip reasons but does not touch ws.nRead/ws.nKept and does not mutate the
  ## record. Callers handle counters, ID synthesis, multiallelic/contig-error
  ## checks (which differ between the hts-nim and C-API paths), and writes.
  result.ok = false

  let svt = resolveSvtype(view, svtypeBuf)
  if svt == svUNKNOWN:
    warnSkip(skUnresolvableSvtype, ws, chrom, pos1, curId,
             "no INFO/SVTYPE and no recognizable ALT")
    return
  if svt == svTRA:
    warnSkip(skUnsupportedTra, ws, chrom, pos1, curId,
             "TRA records are not supported (use BND notation)")
    return
  if svt notin SupportedSvTypes:
    inc ws.skipped[skUnsupportedSvtype]   # INS lands here silently
    return

  if svt == svBND:
    let alts = view.ALT
    let altStr = if alts.len > 0: alts[0] else: ""
    let bnd = parseBndAlt(altStr)
    if not bnd.ok:
      warnSkip(skMalformedBnd, ws, chrom, pos1, curId,
               "could not parse breakend ALT: '" & altStr & "'")
      return
    return (ok: true, svt: svt, endPos: pos1 + 1, svlen: 0'i64, binIdx: 0,
            bndChr2: bnd.chr2, bndPos2: bnd.pos2)

  let endR = resolveEnd(view, endBuf, svlenBuf)
  if not endR.ok:
    warnSkip(skUnresolvableEnd, ws, chrom, pos1, curId,
             "missing both INFO/END and INFO/SVLEN")
    return
  if endR.endPos <= pos1:
    warnSkip(skEndLePos, ws, chrom, pos1, curId,
             "END=" & $endR.endPos & " <= POS=" & $pos1)
    return
  let svR = resolveSvlen(view, endR.endPos, svlenBuf)
  if not svR.ok:
    warnSkip(skUnresolvableEnd, ws, chrom, pos1, curId,
             "could not resolve SVLEN")
    return

  (ok: true, svt: svt, endPos: endR.endPos, svlen: svR.svlen,
   binIdx: binIndexFor(svR.svlen), bndChr2: "", bndPos2: 0'i64)

proc ensureSvInfoDefs(h: vcf.Header) =
  for def in SvInfoDefs:
    let (id, num, typ, desc) = def
    try:
      discard h.get(id, BCF_HEADER_TYPE.BCF_HL_INFO)
    except KeyError:
      discard h.add_info(id, num, typ, desc)

proc normalizeRecord*(hdr: ptr bcf_hdr_t; rec: ptr bcf1_t; lineno: int;
                     ws: var WarnState; view: Variant;
                     svtypeBuf: var string; endBuf, svlenBuf: var seq[int32];
                     vcfPath: string):
    tuple[ok: bool; svt: SvType; endPos: int64; svlen: int64; binIdx: int;
          bndChr2: string; bndPos2: int64] =
  ## Resolve SVTYPE/END/SVLEN/binIdx for a single record. May synthesize an
  ## ID on the record (via bcf_update_id) and may quit() on multiallelic input.
  ## Increments `ws` counters for skipped/synthesized records.
  ## Returns ok=false for any skip reason; caller should not write the record.
  ##
  ## `view` is a reusable non-owning Variant wrapper (from synced_bcf_reader's
  ## newVariantView) that lets us call the Variant-based resolvers on records
  ## owned by a synced reader. Its lifetime must outlive this call.

  inc ws.nRead
  result.ok = false

  if rec.errcode == BCF_ERR_CTG_UNDEF:
    let chrom = getChromName(hdr, rec.rid)
    warnSkip(skUnknownContig, ws, chrom, rec.pos + 1, "", "BCF_ERR_CTG_UNDEF")
    return

  if rec.n_allele.int > 2:
    let chrom = getChromName(hdr, rec.rid)
    logError("multiallelic record at " & chrom & ":" &
             $(rec.pos + 1) & " in " & vcfPath &
             " — split multiallelics before running matcha")
    quit(1)

  discard bcf_unpack(rec, BCF_UN_SHR.cint)

  setRecView(view, hdr, rec)

  let chrom = getChromName(hdr, rec.rid)
  let pos1  = rec.pos + 1
  let curId = getRecId(rec)

  let r = resolveRecord(view, ws, chrom, pos1, curId,
                        svtypeBuf, endBuf, svlenBuf)
  if not r.ok: return

  # Synthesize ID if missing
  if curId.len == 0:
    let newId = synthesizeId(chrom, pos1, r.svt, lineno)
    discard bcf_update_id(hdr, rec, newId.cstring)
    inc ws.syntheticId

  inc ws.nKept
  (ok: true, svt: r.svt, endPos: r.endPos, svlen: r.svlen, binIdx: r.binIdx,
   bndChr2: r.bndChr2, bndPos2: r.bndPos2)

proc emitSummary*(ws: WarnState) =
  let nSkipped = ws.nRead - ws.nKept
  logWarn(ws.callset & " summary: read " & $ws.nRead &
          " records, kept " & $ws.nKept & ", skipped " & $nSkipped)
  var parts: seq[string]
  for r in SkipReason:
    parts.add(reasonStr(r) & "=" & $ws.skipped[r])
  logWarn("  by reason: " & parts.join(", "))
  if ws.syntheticId > 0:
    logWarn("  ids synthesized: " & $ws.syntheticId)

proc tempBcfPath(tmpDir, prefix, svtype: string, binIdx: int): string =
  tmpDir / "matcha_" & $getCurrentProcessId() & "_" & prefix & "_" &
           svtype & "_b" & $binIdx & ".bcf"

proc makeRunTmpDir*(parent: string): string =
  ## Create a unique per-invocation subdirectory under `parent`. Uses mkdtemp
  ## semantics so it is safe against PID reuse across nodes sharing scratch.
  if not dirExists(parent):
    createDir(parent)
  result = createTempDir("matcha_", "", parent)
  logInfo("tmp dir: " & result)

proc hrecToLine(h: ptr bcf_hrec_t): string =
  let keys = cast[ptr UncheckedArray[cstring]](h.keys)
  let vals = cast[ptr UncheckedArray[cstring]](h.vals)
  result = "##INFO=<"
  for i in 0 ..< h.nkeys.int:
    if i > 0: result &= ","
    result &= $keys[i] & "=" & $vals[i]
  result &= ">"

proc buildSlimHdr(src: ptr bcf_hdr_t,
                  infoDefs: openArray[(string, string, string, string)],
                  extraKeepInfo: openArray[string]): ptr bcf_hdr_t =
  ## Duplicate src, strip all FORMAT and INFO defs, then add back only the
  ## fields in infoDefs plus any extraKeepInfo fields (re-serialised verbatim
  ## from the source hrec so that types are preserved for decode in anno).
  result = bcf_hdr_dup(src)
  bcf_hdr_remove(result, BCF_HEADER_TYPE.BCF_HL_FMT.cint, nil)
  bcf_hdr_remove(result, BCF_HEADER_TYPE.BCF_HL_INFO.cint, nil)
  for (id, num, typ, desc) in infoDefs:
    discard bcf_hdr_append(result,
      ("##INFO=<ID=" & id & ",Number=" & num & ",Type=" & typ &
       ",Description=\"" & desc & "\">").cstring)
  if extraKeepInfo.len > 0:
    let extraSet = toHashSet(extraKeepInfo)
    let hrecs = cast[ptr UncheckedArray[ptr bcf_hrec_t]](src.hrec)
    for i in 0 ..< src.nhrec.int:
      let hrec = hrecs[i]
      if hrec.`type` != BCF_HEADER_TYPE.BCF_HL_INFO.cint: continue
      let keys = cast[ptr UncheckedArray[cstring]](hrec.keys)
      let vals = cast[ptr UncheckedArray[cstring]](hrec.vals)
      for j in 0 ..< hrec.nkeys.int:
        if $keys[j] == "ID" and $vals[j] in extraSet:
          discard bcf_hdr_append(result, hrecToLine(hrec).cstring)
          break
  discard bcf_hdr_sync(result)

proc captureChromOrder*(h: vcf.Header): seq[string] =
  ## Return contig IDs in the order they appear in the VCF header.
  var n: cint = 0
  let names = bcf_hdr_seqnames(h.hdr, n.addr)
  if names == nil: return
  for i in 0 ..< n.int:
    result.add($names[i])
  free(names)   # free the array but not the underlying strings (per htslib docs)

proc preprocessVcf*(vcfPath, tmpDir, prefix: string,
                    extraKeepInfo: openArray[string] = [],
                    ioThreads:    int  = 0;
                    noIndex:      bool = false): PreprocOutput =
  ## Stream vcfPath, normalize each record, and write per-SVTYPE BCFs.
  ## When noIndex=false (default), all temp BCFs are CSI-indexed on return.
  ## Inputs may be VCF.gz or BCF.
  ##
  ## extraKeepInfo: additional INFO field names that survive the slim step.
  ## Used by `matcha anno` to carry user-requested DB fields into the
  ## per-(svtype, bin) BCFs alongside the default keep-set.
  logVerbose("[" & prefix & "] reading " & vcfPath)
  var vcf: VCF
  if not open(vcf, vcfPath, threads = ioThreads):
    raise newException(IOError, "cannot open VCF/BCF: " & vcfPath)
  vcf.set_samples(@["^"])
  ensureSvInfoDefs(vcf.header)

  # Per-SVTYPE keep-sets (HashSet membership is O(1) per INFO field per record).
  var keepSetInterval = initHashSet[string]()
  var keepSetBnd      = initHashSet[string]()
  for (id, _, _, _) in IntervalInfoDefs: keepSetInterval.incl(id)
  for (id, _, _, _) in BndInfoDefs:      keepSetBnd.incl(id)
  for n in extraKeepInfo:
    keepSetInterval.incl(n)
    keepSetBnd.incl(n)

  # Slim header templates: built once, duped per writer.
  let slimHdrInterval = buildSlimHdr(vcf.header.hdr, IntervalInfoDefs, extraKeepInfo)
  let slimHdrBnd      = buildSlimHdr(vcf.header.hdr, BndInfoDefs,      extraKeepInfo)

  var writers: Table[SvtypeBin, VCF]
  var indexes: Table[SvtypeBin, ptr hts_idx_t]
  var ws = WarnState(cap: warnCap(), callset: "callset" & prefix)
  var lineno = 0
  var srcIndex: int32 = 0  ## Increments for every record read (including skipped),
                           ## so anno's sequential A-file walk produces matching values.
  var svtypeStr: string
  var endData, svlenData: seq[int32]
  var toDelete: seq[string]

  result.chromOrder = captureChromOrder(vcf.header)

  for v in vcf:
    var curSrcIndex = srcIndex
    inc srcIndex
    inc lineno
    inc ws.nRead

    block recordBody:
      if v.c.errcode == BCF_ERR_CTG_UNDEF:
        warnSkip(skUnknownContig, ws, $v.CHROM, v.POS, $v.ID,
                 "BCF_ERR_CTG_UNDEF")
        break recordBody

      if v.ALT.len > 1:
        logError("multiallelic record at " & $v.CHROM & ":" &
                 $v.POS & " in " & vcfPath &
                 " — split multiallelics before running matcha")
        quit(1)

      let chrom = $v.CHROM
      let pos1  = v.POS
      let curId = $v.ID
      let r = resolveRecord(v, ws, chrom, pos1, curId,
                            svtypeStr, endData, svlenData)
      if not r.ok: break recordBody
      let svt    = r.svt
      let endPos = r.endPos
      let svlen  = r.svlen
      let binIdx = r.binIdx
      let bndChr2 = r.bndChr2
      let bndPos2 = r.bndPos2

      # Synthesize ID if missing
      if curId.len == 0 or curId == ".":
        v.ID = synthesizeId(chrom, pos1, svt, lineno)
        inc ws.syntheticId

      # Slim INFO: drop everything outside the per-SVTYPE keep-set (two-phase
      # to avoid iterator invalidation as we delete).
      let activeKeepSet = if svt == svBND: keepSetBnd else: keepSetInterval
      toDelete.setLen(0)
      for fld in v.info.fields:
        if fld.name notin activeKeepSet:
          toDelete.add(fld.name)
      for name in toDelete:
        discard v.info.delete(name)

      # Authoritative writes — only the fields in the active keep-set.
      if svt == svBND:
        # ALT parse wins over any stale INFO/CHR2 or INFO/POS2 in the input.
        # endPos (= POS+1) is used only for hts_idx_push below; not written.
        var chr2Str = bndChr2
        discard v.info.set("CHR2", chr2Str)
        var pos2I32 = bndPos2.int32
        discard v.info.set("POS2", pos2I32)
      else:
        var endI32 = endPos.int32
        discard v.info.set("END", endI32)

      # Lazily open per-(svtype, bin) writer.
      let key: SvtypeBin = (svt, binIdx)
      if key notin writers:
        let path = tempBcfPath(tmpDir, prefix, $svt, binIdx)
        var wtr: VCF
        if not open(wtr, path, mode = "wb"):
          raise newException(IOError, "cannot create temp BCF: " & path)
        wtr.copy_header(vcf.header)   # initializes wtr.header (open leaves it nil)
        let slimHdr = if svt == svBND: slimHdrBnd else: slimHdrInterval
        bcf_hdr_destroy(wtr.header.hdr)
        wtr.header.hdr = bcf_hdr_dup(slimHdr)
        wtr.set_samples(@["^"])   # temp BCFs carry no sample columns
        discard wtr.write_header()
        writers[key] = wtr
        result.paths[key] = path
        if not noIndex:
          let headerOff = uint64(bgzf_tell(bgzfHandle(wtr)))
          let idx = hts_idx_init(0, HTS_FMT_CSI.cint, headerOff, 14, 5)
          if idx == nil:
            raise newException(IOError, "cannot create CSI index for: " & path)
          indexes[key] = idx
      if svt notin result.chromsBySvtype:
        result.chromsBySvtype[svt] = initHashSet[string]()
      if svt notin result.populatedBins:
        result.populatedBins[svt] = {}
      result.populatedBins[svt].incl(uint8(binIdx))

      # REF/ALT and QUAL are unused by matchcore — blank them to shrink the
      # record.
      discard bcf_update_alleles_str(vcf.header.hdr, v.c, "N\t.")
      v.c.qual = cast[cfloat](bcf_float_missing.uint32)

      # Write the sequential record index (same value whether or not this record
      # was skipped — curSrcIndex matches what anno's counter will see).
      discard v.info.set("SRC_INDEX", curSrcIndex)

      discard writers[key].write_variant(v)
      if not noIndex:
        let woff = uint64(bgzf_tell(bgzfHandle(writers[key])))
        discard hts_idx_push(indexes[key], v.c.rid, int64(v.c.pos), int64(v.c.pos) + int64(v.c.rlen), woff, 1)
      result.chromsBySvtype[svt].incl($v.CHROM)
      inc ws.nKept

  vcf.close()

  for key, wtr in writers.mpairs:
    if not noIndex:
      let finalOff = uint64(bgzf_tell(bgzfHandle(wtr)))
      hts_idx_finish(indexes[key], finalOff)
    wtr.close()
    if not noIndex:
      let path = result.paths[key]
      hts_idx_save(indexes[key], path.cstring, HTS_FMT_CSI.cint)
      hts_idx_destroy(indexes[key])

  bcf_hdr_destroy(slimHdrInterval)
  bcf_hdr_destroy(slimHdrBnd)

  logInfo("[" & prefix & "] wrote " & $result.paths.len &
       " temp BCFs ((svtype, bin) groups)" &
       (if noIndex: " (no index)" else: " (CSI indexed)"))
  emitSummary(ws)

proc buildWorkQueue*(a, b: PreprocOutput,
                     cfg: MatchConfig): tuple[jobs: seq[MatchJob], fileList: seq[string]] =
  ## Emit one MatchJob per (chrom, svtype, binA) where:
  ##   - chrom is in both A and B for this svtype
  ##   - binA exists in A
  ##   - at least one populated B bin is adjacent to binA under the binding
  ##     threshold (max of the supplied --min-overlap / --min-jaccard values).
  ##
  ## Each job carries the path of A's per-(svtype, binA) BCF plus a table
  ## mapping each adjacent populated B bin → (B's path, B's file index).
  ##
  ## Also returns fileList: a deduplicated ordered list of all slim-BCF paths.
  ## FILE_IDX values in MatchJob fields are indices into this list.
  ##
  ## Job order: chrom in input header order, then svtype string, then binA.

  let threshold = cfg.threshold

  # Build the global file list from all A and B paths in sorted key order.
  var fileList: seq[string]
  var pathToIdx: Table[string, int16]
  for key in a.paths.keys.toSeq.sorted:
    let p = a.paths[key]
    if p notin pathToIdx:
      pathToIdx[p] = fileList.len.int16
      fileList.add(p)
  for key in b.paths.keys.toSeq.sorted:
    let p = b.paths[key]
    if p notin pathToIdx:
      pathToIdx[p] = fileList.len.int16
      fileList.add(p)

  # Index of chrom → header order for sorting.
  var chromIdx: Table[string, int]
  for i, c in a.chromOrder:
    chromIdx[c] = i

  var jobs: seq[MatchJob]
  for key, pathA in a.paths:
    let (svt, binA) = key
    if svt notin b.populatedBins: continue

    let adj = adjacentBins(binA, threshold, b.populatedBins[svt])
    if adj.len == 0: continue

    # In self mode, restrict binsB to binB >= binA so each cross-bin pair is
    # processed once. Same-bin pairs are deduped in matchcore via srcIndexA < srcIndexB.
    var binsB: Table[int, BinBEntry]
    for binB in adj:
      if cfg.selfMode and binB < binA: continue
      let bKey: SvtypeBin = (svt, binB)
      if bKey in b.paths:
        let bPath = b.paths[bKey]
        binsB[binB] = (path: bPath, fileIdx: pathToIdx[bPath])
    if binsB.len == 0: continue

    let chromsA = a.chromsBySvtype[svt]
    let chromsB = b.chromsBySvtype[svt]
    for chrom in chromsA:
      if chrom notin chromsB: continue
      jobs.add(MatchJob(
        chrom:    chrom, chromIdx: int16(chromIdx.getOrDefault(chrom, 0)),
        svtype:   svt, binA: binA,
        pathA:    pathA, fileIdxA: pathToIdx[pathA], binsB: binsB,
      ))

  # Stable sort: chrom (header order), then svtype string, then binA.
  proc cmpJobs(x, y: MatchJob): int =
    let xi = chromIdx.getOrDefault(x.chrom, high(int))
    let yi = chromIdx.getOrDefault(y.chrom, high(int))
    if xi != yi: return cmp(xi, yi)
    let s = cmp($x.svtype, $y.svtype)
    if s != 0: return s
    cmp(x.binA, y.binA)

  jobs.sort(cmpJobs)
  result = (jobs: jobs, fileList: fileList)

proc removeTempBcfs*(paths: Table[SvtypeBin, string]) =
  ## Delete per-(svtype, bin) BCF temp files and their CSI indexes.
  for path in paths.values:
    if fileExists(path):          removeFile(path)
    if fileExists(path & ".csi"): removeFile(path & ".csi")

# ---------------------------------------------------------------------------
# Parallel preprocessing helper (shared by match and anno)
# ---------------------------------------------------------------------------

type PreprocInput* = object
  path*:      string
  tmpDir*:    string
  prefix*:    string
  extraKeep*: seq[string]
  ioThreads*: int

var gPpIn:  array[2, PreprocInput]
var gPpOut: array[2, PreprocOutput]

proc ppWorker(idx: int) {.thread.} =
  {.cast(gcsafe).}:
    let s = gPpIn[idx]
    gPpOut[idx] = preprocessVcf(s.path, s.tmpDir, s.prefix, s.extraKeep, s.ioThreads)

proc runParallelPreproc*(a, b: PreprocInput): tuple[a, b: PreprocOutput] =
  gPpIn[0] = a; gPpIn[1] = b
  var thA, thB: Thread[int]
  createThread(thA, ppWorker, 0)
  createThread(thB, ppWorker, 1)
  joinThread(thA); joinThread(thB)
  (gPpOut[0], gPpOut[1])
