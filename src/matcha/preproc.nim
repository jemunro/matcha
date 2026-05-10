## preproc.nim — VCF/BCF normalization, slimming, and per-(SVTYPE, bin) BCF output.
##
## For each input record we:
##   * Resolve SVTYPE (INFO/SVTYPE; symbolic ALT fallback; ALT wins on conflict).
##   * Resolve END   (INFO/END; POS+abs(SVLEN) fallback; require END > POS).
##   * Resolve SVLEN (INFO/SVLEN abs; END-POS fallback; warn+normalize on >10%
##     disagreement when both INFO fields were independently provided).
##   * Synthesize ID if absent ("." or empty) → CHROM_POS_SVTYPE_LINENUMBER.
##   * Slim INFO to the keep-set (SVTYPE, SVLEN, END, CHR2, END2, POS2,
##     MATCHA_BOFF).
##   * Compute size-bin via bins.binIndexFor(SVLEN).
##   * Write to a per-(SVTYPE, bin) BCF spanning all chroms; CSI-indexed.
##
## Skip categories (counted; final summary always emitted):
##   skUnsupportedSvtype  - BND/INS/TRA (count only, awaiting later milestone)
##   skUnresolvableSvtype - neither INFO/SVTYPE nor symbolic ALT
##   skUnresolvableEnd    - missing both INFO/END and INFO/SVLEN
##   skEndLePos           - resolved END <= POS
##   skUnknownContig      - BCF_ERR_CTG_UNDEF on read
##
## Per-record warnings are throttled at MATCHA_WARN_CAP (default 5) per reason.

import std/[algorithm, os, sets, strutils, tables]
import hts
import hts/private/hts_concat  # BCF_ERR_CTG_UNDEF, BGZF, htsFile, bgzf_tell
import utils, log, bins

# hts-nim's VCF type keeps the underlying htsFile* private. We need access to
# it (specifically vcf.hts.fp.bgzf) so we can call bgzf_tell to record each
# record's source-file BGZF virtual offset during preprocessing. Rather than
# patching vendored hts-nim, mirror enough of VCF's prefix to cast across.
# Both inherit RootObj, so the type-info slot lines up; the very next slot is
# `hts: ptr htsFile`. Layout is documented at vendor/hts-nim/src/hts/vcf.nim:19-26.
type
  VcfPriv = ref object of RootObj
    hts: ptr htsFile

proc bgzfHandle(v: VCF): ptr BGZF {.inline.} =
  cast[VcfPriv](v).hts.fp.bgzf

type
  SvtypeBin* = tuple[svtype: SvType, bin: int]

  PreprocOutput* = object
    paths*:          Table[SvtypeBin, string]         ## (svtype, bin) → BCF path
    populatedBins*:  Table[SvType, set[uint8]]        ## svtype → set of bin indexes
    chromsBySvtype*: Table[SvType, HashSet[string]]   ## svtype → chroms seen
    chromOrder*:     seq[string]                      ## chroms in input header order

  MatchJob* = object
    chrom*:   string
    svtype*:  SvType
    binA*:    int
    pathA*:   string
    binsB*:   Table[int, string]                      ## binB → B's path

  SkipReason = enum
    skUnsupportedSvtype
    skUnresolvableSvtype
    skUnresolvableEnd
    skEndLePos
    skUnknownContig

  WarnState = object
    skipped:      array[SkipReason, int]
    emitted:      array[SkipReason, int]
    inconsistent: int
    syntheticId:  int
    nRead:        int
    nKept:        int
    cap:          int
    callset:      string

const KeepInfo = ["SVTYPE", "SVLEN", "END", "CHR2", "END2", "POS2", "MATCHA_BOFF"]

# Standard SV INFO definitions to backfill when the input header omits them.
# Without these, info.set() at write time fails because hts-nim resolves
# the field's type via the (reader's) header.
const SvInfoDefs = [
  ("SVTYPE",      "1", "String",  "Type of structural variant"),
  ("SVLEN",       "1", "Integer", "Length of the SV (absolute value, positive)"),
  ("END",         "1", "Integer", "End position of the SV (1-based, inclusive)"),
  ("CHR2",        "1", "String",  "Chromosome of mate breakend (BND/TRA)"),
  ("POS2",        "1", "Integer", "Position of mate breakend (BND/TRA)"),
  ("MATCHA_BOFF", "2", "Integer",
   "matcha-internal: source-file BGZF virtual offset (high32, low32)"),
]

proc reasonStr(r: SkipReason): string =
  case r
  of skUnsupportedSvtype:   "unsupported_svtype"
  of skUnresolvableSvtype:  "unresolvable_svtype"
  of skUnresolvableEnd:     "unresolvable_end"
  of skEndLePos:            "end_le_pos"
  of skUnknownContig:       "unknown_contig"

proc warnSkip(reason: SkipReason, ws: var WarnState,
              chrom: string, pos: int64, id: string, detail: string) =
  inc ws.skipped[reason]
  if ws.emitted[reason] < ws.cap:
    inc ws.emitted[reason]
    logWarn(chrom & ":" & $pos & " ID=" & id &
            " reason=\"" & reasonStr(reason) & ": " & detail & "\"")

proc warnInconsistency(ws: var WarnState,
                       chrom: string, pos: int64, id: string, detail: string) =
  inc ws.inconsistent
  if ws.inconsistent <= ws.cap:
    logWarn(chrom & ":" & $pos & " ID=" & id &
            " reason=\"end_svlen_inconsistent: " & detail & "\"")

proc symbolicAltSvtype(alt: string): SvType =
  ## Parse `<DEL>`, `<DUP>`, `<INV>`, `<INS>`, `<TRA>` → corresponding SvType.
  ## Returns svUNKNOWN for non-symbolic ALTs (e.g. BND breakend strings).
  if alt.len < 3 or alt[0] != '<' or alt[^1] != '>':
    return svUNKNOWN
  parseSvType(alt[1 .. ^2])

proc resolveSvtype(v: Variant, infoBuf: var string): SvType =
  ## Resolve SVTYPE. Symbolic ALT (when present and parseable as a known
  ## SvType) takes priority; INFO/SVTYPE is the fallback. Returns svUNKNOWN
  ## if neither resolves. The "ALT wins" precedence matches the contract
  ## documented in CLAUDE.md.
  let infoOk = v.info().get("SVTYPE", infoBuf) == Status.OK and infoBuf.len > 0
  let infoSv =
    if infoOk: parseSvType($cast[cstring](addr infoBuf[0]))
    else:      svUNKNOWN
  var altSv = svUNKNOWN
  let alts = v.ALT
  if alts.len > 0:
    altSv = symbolicAltSvtype(alts[0])
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

proc inconsistencyExceedsTenPct(svlen, endDerived: int64): bool =
  if svlen <= 0: return false
  abs(svlen - endDerived).float / svlen.float > 0.10

proc synthesizeId(chrom: string, pos: int64, svt: SvType, lineno: int): string =
  chrom & "_" & $pos & "_" & $svt & "_" & $lineno

proc ensureSvInfoDefs(h: vcf.Header) =
  for def in SvInfoDefs:
    let (id, num, typ, desc) = def
    try:
      discard h.get(id, BCF_HEADER_TYPE.BCF_HL_INFO)
    except KeyError:
      discard h.add_info(id, num, typ, desc)

proc emitSummary(ws: WarnState) =
  let nSkipped = ws.nRead - ws.nKept
  logWarn(ws.callset & " summary: read " & $ws.nRead &
          " records, kept " & $ws.nKept & ", skipped " & $nSkipped)
  var parts: seq[string]
  for r in SkipReason:
    parts.add(reasonStr(r) & "=" & $ws.skipped[r])
  logWarn("  by reason: " & parts.join(", "))
  logWarn("  inconsistencies: " & $ws.inconsistent &
          " END/SVLEN >10% (END used as authoritative)")
  logWarn("  ids synthesized: " & $ws.syntheticId)

proc tempBcfPath(tmpDir, prefix, svtype: string, binIdx: int): string =
  tmpDir / "matcha_" & $getCurrentProcessId() & "_" & prefix & "_" &
           svtype & "_b" & $binIdx & ".bcf"

proc captureChromOrder(h: vcf.Header): seq[string] =
  ## Return contig IDs in the order they appear in the VCF header.
  var n: cint = 0
  let names = bcf_hdr_seqnames(h.hdr, n.addr)
  if names == nil: return
  for i in 0 ..< n.int:
    result.add($names[i])
  free(names)   # free the array but not the underlying strings (per htslib docs)

proc preprocessVcf*(vcfPath, tmpDir, prefix: string): PreprocOutput =
  ## Stream vcfPath, normalize each record, and write per-SVTYPE BCFs.
  ## All temp BCFs are CSI-indexed on return. Inputs may be VCF.gz or BCF.
  logV("[" & prefix & "] reading " & vcfPath)
  var vcf: VCF
  if not open(vcf, vcfPath):
    raise newException(IOError, "cannot open VCF/BCF: " & vcfPath)
  vcf.set_samples(@["^"])
  ensureSvInfoDefs(vcf.header)

  var writers: Table[SvtypeBin, VCF]
  var ws = WarnState(cap: warnCap(), callset: "callset" & prefix)
  var lineno = 0
  var svtypeStr: string
  var endData, svlenData: seq[int32]

  result.chromOrder = captureChromOrder(vcf.header)

  # Capture the source-file BGZF virtual offset for each record.
  # `bgzf_tell` at the bottom of the loop body returns the position of the
  # *next* record (because the iterator's bcf_read for record N has already
  # advanced past it before yielding). So we prime with one tell before the
  # loop, then update unconditionally at the bottom. The `block recordBody`
  # makes every skip path flow through that final update.
  var nextOffset = int64(bgzf_tell(bgzfHandle(vcf)))

  for v in vcf:
    let recordOffset = nextOffset
    inc lineno
    inc ws.nRead

    block recordBody:
      if v.c.errcode == BCF_ERR_CTG_UNDEF:
        warnSkip(skUnknownContig, ws, $v.CHROM, v.POS, $v.ID,
                 "BCF_ERR_CTG_UNDEF")
        break recordBody

      let svt = resolveSvtype(v, svtypeStr)
      if svt == svUNKNOWN:
        warnSkip(skUnresolvableSvtype, ws, $v.CHROM, v.POS, $v.ID,
                 "no INFO/SVTYPE and no symbolic ALT")
        break recordBody
      if svt notin SupportedSvTypes:
        inc ws.skipped[skUnsupportedSvtype]
        break recordBody

      let (eOk, endPos, endFromInfo) = resolveEnd(v, endData, svlenData)
      if not eOk:
        warnSkip(skUnresolvableEnd, ws, $v.CHROM, v.POS, $v.ID,
                 "missing both INFO/END and INFO/SVLEN")
        break recordBody

      if endPos <= v.POS:
        warnSkip(skEndLePos, ws, $v.CHROM, v.POS, $v.ID,
                 "END=" & $endPos & " <= POS=" & $v.POS)
        break recordBody

      let (sOk, svlenInit, svlenFromInfo) = resolveSvlen(v, endPos, svlenData)
      if not sOk:
        warnSkip(skUnresolvableEnd, ws, $v.CHROM, v.POS, $v.ID,
                 "could not resolve SVLEN")
        break recordBody

      let endDerived = endPos - v.POS
      var svlen = svlenInit
      if endFromInfo and svlenFromInfo and
         inconsistencyExceedsTenPct(svlen, endDerived):
        warnInconsistency(ws, $v.CHROM, v.POS, $v.ID,
                          "INFO/SVLEN=" & $svlen & " vs END-POS=" & $endDerived)
        svlen = endDerived

      # Synthesize ID if missing
      let curId = $v.ID
      if curId.len == 0 or curId == ".":
        v.ID = synthesizeId($v.CHROM, v.POS, svt, lineno)
        inc ws.syntheticId

      # Slim INFO: drop everything outside the keep-set (two-phase to avoid
      # iterator invalidation as we delete).
      var toDelete: seq[string]
      for fld in v.info.fields:
        if fld.name notin KeepInfo:
          toDelete.add(fld.name)
      for name in toDelete:
        discard v.info.delete(name)

      # Authoritative writes (overwrite whatever was there with normalized values)
      var svtStr = $svt
      discard v.info.set("SVTYPE", svtStr)
      var svlenI32 = svlen.int32
      discard v.info.set("SVLEN", svlenI32)
      var endI32 = endPos.int32
      discard v.info.set("END", endI32)

      # Encode the source-file offset as INFO/MATCHA_BOFF (Number=2, Integer)
      var boffPair: seq[int32] = @[
        int32((recordOffset shr 32) and 0xFFFFFFFF'i64),
        int32(recordOffset and 0xFFFFFFFF'i64),
      ]
      discard v.info.set("MATCHA_BOFF", boffPair)

      # Compute bin and lazily open per-(svtype, bin) writer.
      let binIdx = binIndexFor(svlen)
      let key: SvtypeBin = (svt, binIdx)
      if key notin writers:
        let path = tempBcfPath(tmpDir, prefix, $svt, binIdx)
        var wtr: VCF
        if not open(wtr, path, mode = "wb"):
          raise newException(IOError, "cannot create temp BCF: " & path)
        wtr.copy_header(vcf.header)
        discard wtr.write_header()
        writers[key] = wtr
        result.paths[key] = path
      if svt notin result.chromsBySvtype:
        result.chromsBySvtype[svt] = initHashSet[string]()
      if svt notin result.populatedBins:
        result.populatedBins[svt] = {}
      result.populatedBins[svt].incl(uint8(binIdx))

      discard writers[key].write_variant(v)
      result.chromsBySvtype[svt].incl($v.CHROM)
      inc ws.nKept

    # Position is now at the start of the next record; capture it for the
    # next iteration. Runs unconditionally — every continue path inside the
    # body uses `break recordBody` so we always reach this update.
    nextOffset = int64(bgzf_tell(bgzfHandle(vcf)))

  vcf.close()

  for key, wtr in writers.mpairs:
    wtr.close()
    let path = result.paths[key]
    bcfBuildIndex(path, path & ".csi", csi = true, threads = 1)

  logV("[" & prefix & "] indexed " & $result.paths.len &
       " temp BCFs ((svtype, bin) groups)")
  emitSummary(ws)

proc buildWorkQueue*(a, b: PreprocOutput, cfg: MatchConfig): seq[MatchJob] =
  ## Emit one MatchJob per (chrom, svtype, binA) where:
  ##   - chrom is in both A and B for this svtype
  ##   - binA exists in A
  ##   - at least one populated B bin is adjacent to binA under the binding
  ##     threshold (max of the supplied --min-overlap / --min-jaccard values).
  ##
  ## Each job carries the path of A's per-(svtype, binA) BCF plus a table
  ## mapping each adjacent populated B bin → B's per-(svtype, binB) BCF path.
  ##
  ## Job order: chrom in input header order, then svtype string, then binA.

  # Binding threshold for size-bin adjacency: the stricter of the two.
  # (Both metrics share the [lenA*t, lenA/t] size constraint, and a higher
  # threshold yields a smaller adjacent set.)
  let threshold =
    if cfg.minOverlapSet and cfg.minJaccardSet:
      max(cfg.minOverlap, cfg.minJaccard)
    elif cfg.minOverlapSet:
      cfg.minOverlap
    else:
      cfg.minJaccard

  # Index of chrom → header order for sorting. Use A's chromOrder; fall back
  # to alphabetical for any chrom not in A's header (shouldn't happen since
  # A is always the source we group by).
  var chromIdx: Table[string, int]
  for i, c in a.chromOrder:
    chromIdx[c] = i

  var jobs: seq[MatchJob]
  for key, pathA in a.paths:
    let (svt, binA) = key
    if svt notin b.populatedBins: continue

    let adj = adjacentBins(binA, threshold, b.populatedBins[svt])
    if adj.len == 0: continue

    # In self mode, A and B are the same callset. Restrict binsB to
    # binB >= binA so each cross-bin pair (binA, binB) is processed once
    # rather than twice — same-bin pairs still need result-level dedup
    # via aOff < bOff, but cross-bin pairs are naturally one-way here.
    var binsB: Table[int, string]
    for binB in adj:
      if cfg.selfMode and binB < binA: continue
      let bKey: SvtypeBin = (svt, binB)
      if bKey in b.paths:
        binsB[binB] = b.paths[bKey]
    if binsB.len == 0: continue

    let chromsA = a.chromsBySvtype[svt]
    let chromsB = b.chromsBySvtype[svt]
    for chrom in chromsA:
      if chrom notin chromsB: continue
      jobs.add(MatchJob(
        chrom:  chrom, svtype: svt, binA: binA,
        pathA:  pathA, binsB: binsB,
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
  result = jobs
