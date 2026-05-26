## matchcore.nim — shared per-job matching loop.
##
## Both `matcha match`, `matcha anno`, and `matcha collapse` walk the same
## A-vs-B candidate space. This module is intentionally minimal: it returns
## MatchPair triples and leaves per-mode field resolution to the adapters.
##
## Self-mode dedup (srcIndexA < srcIndexB filter) is baked in via cfg.selfMode.
## Singleton emission (A records with no passing B match) is controlled by
## cfg.emitSingletons; only collapse sets this true.

import std/[algorithm, deques, math, sequtils, tables]
import std/options as stdopt
import hts
import hts/private/hts_concat
import intervals, preproc, bins, utils

# ---------------------------------------------------------------------------
# tearDownVcf — close a VCF and free the htslib header eagerly.
#
# hts-nim's `close()` only closes the file handle; the BCF header, CSI/tabix
# indexes, and `bcf1_t` buffer are left for the GC finalizer to free. Under
# ORC the finalizer is not guaranteed to run between jobs, so opening many
# per-bin BCFs per match job accumulates C-heap (header is by far the largest
# piece — ~500 KB per VCF with N samples × INFO/FORMAT lines).
#
# Freeing `header.hdr` here drops the per-job leak ~50×. `bidx`, `tidx`, and
# `c` are private fields of hts-nim's VCF and still leak ~0.2 MB per opened
# VCF until the GC eventually runs; bounded and small enough to accept rather
# than reach into private layout.
# ---------------------------------------------------------------------------

proc tearDownVcf*(v: VCF) =
  if v == nil: return
  v.close()
  if v.header != nil and v.header.hdr != nil:
    bcf_hdr_destroy(v.header.hdr)
    v.header.hdr = nil

# ---------------------------------------------------------------------------
# Slim-BCF INFO decode helpers
# ---------------------------------------------------------------------------

proc passFilterId*(v: VCF): cint {.inline.} =
  ## FILTER dictionary ID for "PASS" in this VCF's header. Cache once per VCF
  ## handle to avoid per-record string allocation by hts-nim's FILTER getter.
  bcf_hdr_id2int(v.header.hdr, BCF_DT_ID.cint, "PASS".cstring)

proc isPass*(va: Variant; passId: cint): bool {.inline.} =
  ## True iff FILTER is unset (n_flt == 0, treated as PASS by hts-nim) or
  ## a single entry equal to PASS. Matches `$va.FILTER == "PASS"` exactly
  ## without the string allocation.
  let n = va.c.d.n_flt
  if n == 0: return true
  if n == 1:
    let arr = cast[ptr UncheckedArray[cint]](va.c.d.flt)
    return arr[0] == passId
  false

proc readSrcIndex*(v: Variant; scratch: var seq[int32]): int32 =
  ## Decode INFO/SRC_INDEX (Number=1 Integer). Returns -1 if absent.
  if v.info().get("SRC_INDEX", scratch) != Status.OK or scratch.len < 1:
    return -1'i32
  scratch[0]

proc readPos2*(v: Variant; scratch: var seq[int32]): tuple[ok: bool; pos2: int64] =
  ## Decode INFO/POS2 (Number=1 Integer). Returns ok=false if absent.
  if v.info().get("POS2", scratch) != Status.OK or scratch.len < 1:
    return (false, 0'i64)
  (true, int64(scratch[0]))

proc readChr2*(v: Variant; scratch: var string): tuple[ok: bool; chr2: string] =
  ## Decode INFO/CHR2 (Number=1 String). Returns ok=false if absent or empty.
  if v.info().get("CHR2", scratch) != Status.OK or scratch.len == 0:
    return (false, "")
  (true, scratch)

proc readSvlen*(v: Variant; scratch: var seq[int32]): tuple[ok: bool; svlen: int64] =
  ## Decode INFO/SVLEN (Number=1 Integer, abs). Returns ok=false if absent or 0.
  if v.info().get("SVLEN", scratch) != Status.OK or scratch.len < 1:
    return (false, 0'i64)
  let n = abs(int64(scratch[0]))
  if n <= 0: return (false, 0'i64)
  (true, n)

proc extractEnd*(v: Variant; endData, svlenData: var seq[int32];
                 outEnd: var int64): bool =
  ## Resolve END for a slim record. After preproc, END is always written
  ## authoritatively, so this is a fast path; SVLEN fallback is defensive.
  if v.info().get("END", endData) == Status.OK and endData.len > 0:
    outEnd = int64(endData[0]); return true
  if v.info().get("SVLEN", svlenData) == Status.OK and svlenData.len > 0:
    outEnd = v.POS + abs(int64(svlenData[0])); return true
  false

# ---------------------------------------------------------------------------
# Sliding-cache helper — shared by BND and INS streamers
# ---------------------------------------------------------------------------

proc advanceSlidingCache*[T](
    vcfB: VCF; chrom: string; slop: int; posA: int64;
    cache: var Deque[T]; cacheEnd: var int64;
    decodeB: proc(vb: Variant): stdopt.Option[T] {.closure.}) =
  ## Maintain a position-ordered sliding cache of B records covering the
  ## window [posA - slop + 1, posA + slop). Evicts the front below the
  ## window and lazily fetches the right edge via a single CSI region query.
  ## `T` must expose `pos: int64`; `decodeB` returns `none` to skip a record.
  let winLo = posA - slop.int64 + 1
  let winHi = posA + slop.int64
  while cache.len > 0 and cache.peekFirst().pos < winLo:
    cache.popFirst()
  let queryLo = max(cacheEnd, winLo)
  if queryLo < winHi:
    let regLo = max(1'i64, queryLo)
    let regHi = winHi - 1
    if regHi >= regLo:
      let region = chrom & ":" & $regLo & "-" & $regHi
      for vb in vcfB.query(region):
        if vb.POS < queryLo: continue
        let rec = decodeB(vb)
        if stdopt.isSome(rec): cache.addLast(stdopt.get(rec))
    cacheEnd = winHi

# ---------------------------------------------------------------------------
# Interval matching
# ---------------------------------------------------------------------------

proc singletonPair(job: MatchJob; srcIndex: int32; posA: int64;
                   passA: bool; qualQ: uint16; callerIdxA: int16): MatchPair =
  MatchPair(srcIndexA: srcIndex, srcIndexB: NO_MATCH,
            posA: int32(posA), posB: 0, sim: 0.0f32,
            fileIdxA: job.fileIdxA, fileIdxB: int16(NO_MATCH),
            chromIdx: job.chromIdx, svtype: int8(job.svtype),
            passA: passA, qualQ: qualQ, callerIdxA: callerIdxA)

proc streamJobPairs*(job: MatchJob; cfg: MatchConfig): seq[MatchPair] =
  ## Drive the per-job match: stream A from the per-(svtype, binA) BCF
  ## restricted to job.chrom, fetch B candidates lazily through TiledBuffers,
  ## compute the active interval metric, apply the threshold, and append a
  ## MatchPair for each passing pair.
  ## When cfg.emitSingletons, also emit a singleton MatchPair for every A
  ## record that has no passing B match.
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

  var endData, svlenData, idxData, ciData: seq[int32]
  let passId = passFilterId(vcfA)
  let metricIsOverlap = cfg.metric == mOverlap

  proc fetchTile(binB, tileIdx: int): seq[BufferedRec] =
    if binB notin vcfsB:
      var v: VCF
      if not open(v, job.binsB[binB].path):
        raise newException(IOError, "cannot open B BCF: " & job.binsB[binB].path)
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
        srcIndex: readSrcIndex(vb, idxData),
        fileIdx:  job.binsB[binB].fileIdx,
      ))

  for va in vcfA.query(job.chrom & ":" & $job.posStart & "-" & $job.posEnd):
    var endA: int64
    if not extractEnd(va, endData, svlenData, endA): continue
    let posA       = va.POS
    let srcIndexA  = readSrcIndex(va, idxData)
    let passA      = isPass(va, passId)
    let qualQ      = quantizeQual(va.QUAL.float32)
    let callerIdxA = (if va.info().get("CALLER_IDX", ciData) == Status.OK and
                        ciData.len > 0: int16(ciData[0]) else: 0'i16)

    var anyMatch = false
    for binB in sortedBinsB:
      let b = binB
      let cands = buffers[b].getCandidates(posA, endA,
        proc(ti: int): seq[BufferedRec] = fetchTile(b, ti))
      for cand in cands:
        if cfg.selfMode and srcIndexA >= cand.srcIndex: continue
        let sim =
          if metricIsOverlap:
            reciprocalOverlap(posA, endA, cand.pos, cand.endPos)
          else:
            jaccard(posA, endA, cand.pos, cand.endPos)
        if sim < cfg.threshold: continue
        result.add(MatchPair(
          srcIndexA: srcIndexA, srcIndexB: cand.srcIndex,
          posA:      int32(posA), posB: int32(cand.pos),
          sim:       float32(sim),
          fileIdxA:  job.fileIdxA, fileIdxB: cand.fileIdx,
          chromIdx:  job.chromIdx, svtype: int8(job.svtype),
          passA: passA, qualQ: qualQ, callerIdxA: callerIdxA,
        ))
        anyMatch = true

    if cfg.emitSingletons and not anyMatch:
      result.add(singletonPair(job, srcIndexA, posA, passA, qualQ, callerIdxA))

    for binB, buf in buffers.mpairs:
      buf.evict(posA)

  tearDownVcf(vcfA)
  for v in vcfsB.mvalues:
    tearDownVcf(v)

# ---------------------------------------------------------------------------
# BND matching (point events; slop-based proximity)
# ---------------------------------------------------------------------------

proc streamBndJobPairs*(job: MatchJob; cfg: MatchConfig): seq[MatchPair] =
  ## Per-job BND matching with a sliding cache of B records.
  let slop = cfg.bndSlop
  if slop <= 0: return
  if 0 notin job.binsB: return

  var vcfA: VCF
  if not open(vcfA, job.pathA):
    raise newException(IOError, "cannot open A BND BCF: " & job.pathA)
  var vcfB: VCF
  if not open(vcfB, job.binsB[0].path):
    raise newException(IOError, "cannot open B BND BCF: " & job.binsB[0].path)

  type BndCacheRec = object
    pos:      int64
    pos2:     int64
    srcIndex: int32
    fileIdx:  int16
    chr2:     string

  var cache = initDeque[BndCacheRec]()
  var cacheEnd: int64 = low(int64)

  var idxData, pos2Data, ciData: seq[int32]
  var chr2Data: string
  let twoSlop = float64(2 * slop)
  let bFileIdx = job.binsB[0].fileIdx
  let passId = passFilterId(vcfA)

  let decodeB = proc(vb: Variant): stdopt.Option[BndCacheRec] =
    let p2B = readPos2(vb, pos2Data)
    if not p2B.ok: return stdopt.none(BndCacheRec)
    let c2B = readChr2(vb, chr2Data)
    if not c2B.ok: return stdopt.none(BndCacheRec)
    stdopt.some(BndCacheRec(
      pos: vb.POS, pos2: p2B.pos2,
      srcIndex: readSrcIndex(vb, idxData),
      fileIdx:  bFileIdx,
      chr2:     c2B.chr2))

  for va in vcfA.query(job.chrom & ":" & $job.posStart & "-" & $job.posEnd):
    let posA       = va.POS
    let srcIndexA  = readSrcIndex(va, idxData)
    let p2A = readPos2(va, pos2Data)
    if not p2A.ok: continue
    let c2A = readChr2(va, chr2Data)
    if not c2A.ok: continue
    let pos2A = p2A.pos2
    let chr2A = c2A.chr2
    let passA      = isPass(va, passId)
    let qualQ      = quantizeQual(va.QUAL.float32)
    let callerIdxA = (if va.info().get("CALLER_IDX", ciData) == Status.OK and
                        ciData.len > 0: int16(ciData[0]) else: 0'i16)

    advanceSlidingCache(vcfB, job.chrom, slop, posA, cache, cacheEnd, decodeB)

    var anyMatch = false
    for i in 0 ..< cache.len:
      let cand = cache[i]
      if cand.chr2 != chr2A: continue
      let d2 = abs(cand.pos2 - pos2A)
      if d2 >= slop: continue
      let d1 = abs(cand.pos - posA)
      let sim = (twoSlop - float64(d1) - float64(d2)) / twoSlop
      if sim <= 0: continue
      if cfg.selfMode and srcIndexA >= cand.srcIndex: continue
      result.add(MatchPair(
        srcIndexA: srcIndexA, srcIndexB: cand.srcIndex,
        posA:      int32(posA), posB: int32(cand.pos),
        sim:       float32(sim),
        fileIdxA:  job.fileIdxA, fileIdxB: cand.fileIdx,
        chromIdx:  job.chromIdx, svtype: int8(job.svtype),
        passA: passA, qualQ: qualQ, callerIdxA: callerIdxA,
      ))
      anyMatch = true

    if cfg.emitSingletons and not anyMatch:
      result.add(singletonPair(job, srcIndexA, posA, passA, qualQ, callerIdxA))

  tearDownVcf(vcfA)
  tearDownVcf(vcfB)

# ---------------------------------------------------------------------------
# INS matching (point events; combined position + size similarity)
# ---------------------------------------------------------------------------

proc streamInsJobPairs*(job: MatchJob; cfg: MatchConfig): seq[MatchPair] =
  ## Per-job INS matching with a sliding cache of B records, mirroring the
  ## BND streamer. INS records all share bin 0; combined similarity is
  ##   sim = sqrt(pos_sim * len_sim)
  ## where pos_sim = (slop - |dPos|) / slop and len_sim = min/max(SVLEN).
  let slop = cfg.insSlop
  let minSim = cfg.insMinSim
  if slop <= 0: return
  if 0 notin job.binsB: return

  var vcfA: VCF
  if not open(vcfA, job.pathA):
    raise newException(IOError, "cannot open A INS BCF: " & job.pathA)
  var vcfB: VCF
  if not open(vcfB, job.binsB[0].path):
    raise newException(IOError, "cannot open B INS BCF: " & job.binsB[0].path)

  type InsCacheRec = object
    pos:      int64
    svlen:    int64
    srcIndex: int32
    fileIdx:  int16

  var cache = initDeque[InsCacheRec]()
  var cacheEnd: int64 = low(int64)

  var idxData, svlenData, ciData: seq[int32]
  let slopF = slop.float64
  let bFileIdx = job.binsB[0].fileIdx
  let passId = passFilterId(vcfA)

  let decodeB = proc(vb: Variant): stdopt.Option[InsCacheRec] =
    let svB = readSvlen(vb, svlenData)
    if not svB.ok: return stdopt.none(InsCacheRec)
    stdopt.some(InsCacheRec(
      pos:      vb.POS,
      svlen:    svB.svlen,
      srcIndex: readSrcIndex(vb, idxData),
      fileIdx:  bFileIdx))

  for va in vcfA.query(job.chrom & ":" & $job.posStart & "-" & $job.posEnd):
    let posA       = va.POS
    let srcIndexA  = readSrcIndex(va, idxData)
    let svA = readSvlen(va, svlenData)
    if not svA.ok: continue
    let svlenA = svA.svlen
    let passA      = isPass(va, passId)
    let qualQ      = quantizeQual(va.QUAL.float32)
    let callerIdxA = (if va.info().get("CALLER_IDX", ciData) == Status.OK and
                        ciData.len > 0: int16(ciData[0]) else: 0'i16)

    advanceSlidingCache(vcfB, job.chrom, slop, posA, cache, cacheEnd, decodeB)

    var anyMatch = false
    for i in 0 ..< cache.len:
      let cand = cache[i]
      let dPos = abs(cand.pos - posA)
      if dPos >= slop: continue
      let posSim = (slopF - dPos.float64) / slopF
      let svF1 = svlenA.float64
      let svF2 = cand.svlen.float64
      let lenSim = min(svF1, svF2) / max(svF1, svF2)
      let sim = sqrt(posSim * lenSim)
      if sim < minSim: continue
      if cfg.selfMode and srcIndexA >= cand.srcIndex: continue
      result.add(MatchPair(
        srcIndexA: srcIndexA, srcIndexB: cand.srcIndex,
        posA:      int32(posA), posB: int32(cand.pos),
        sim:       float32(sim),
        fileIdxA:  job.fileIdxA, fileIdxB: cand.fileIdx,
        chromIdx:  job.chromIdx, svtype: int8(job.svtype),
        passA: passA, qualQ: qualQ, callerIdxA: callerIdxA,
      ))
      anyMatch = true

    if cfg.emitSingletons and not anyMatch:
      result.add(singletonPair(job, srcIndexA, posA, passA, qualQ, callerIdxA))

  tearDownVcf(vcfA)
  tearDownVcf(vcfB)
