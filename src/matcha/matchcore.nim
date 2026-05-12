## matchcore.nim — shared per-job matching loop.
##
## Both `matcha match` and `matcha anno` walk the same A-vs-B candidate space.
## The differences are narrow:
##   - what (if anything) gets pulled off each B candidate, and
##   - what gets emitted (and whether) for each passing pair.
##
## streamJobPairs[B, R] handles the generic streaming + tiled-buffer + filter
## machinery and invokes the supplied `extract` and `emit` callbacks at the
## right points. Per-mode adapters live in match.nim and anno.nim.

import std/[algorithm, deques, sequtils, tables]
import hts
import intervals, preproc, bins, utils

# ---------------------------------------------------------------------------
# Shared low-level helpers (previously duplicated in match.nim)
# ---------------------------------------------------------------------------

proc readBoff*(v: Variant; scratch: var seq[int32]): int64 =
  ## Decode INFO/MATCHA_BOFF (Number=2 Integer: high32, low32) into an int64.
  ## Returns 0 if the field is absent. The mask on the low half avoids sign-
  ## extension when the int32 has its high bit set.
  if v.info().get("MATCHA_BOFF", scratch) != Status.OK or scratch.len < 2:
    return 0
  result = (int64(scratch[0]) shl 32) or (int64(scratch[1]) and 0xFFFFFFFF'i64)

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
# Generic per-job streaming
# ---------------------------------------------------------------------------

type
  PairResult*[R] = object
    keep*:  bool   ## false → callback rejected this pair (e.g. self-mode dedup)
    item*:  R

  BExtractCb*[B] = proc(v: Variant): B {.closure.}
  PairEmitCb*[B, R] = proc(va: Variant; posA, endA, aOff: int64;
                            cand: BufferedRec; bExtra: B;
                            sim: float64): PairResult[R] {.closure.}

proc streamJobPairs*[B, R](job: MatchJob; cfg: MatchConfig;
                           extract: BExtractCb[B];
                           emit: PairEmitCb[B, R]): seq[R] =
  ## Drive the per-job match: stream A from the per-(svtype, binA) BCF
  ## restricted to job.chrom, fetch B candidates lazily through TiledBuffers,
  ## compute the active interval metric (reciprocal overlap or Jaccard,
  ## chosen via cfg.metric), apply the threshold, and dispatch each passing
  ## pair to `emit`. The optional `extract` callback runs once per fetched
  ## B record and its return value is threaded into the matching `emit`
  ## call — used by anno mode to carry user-requested INFO values.
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

  # Parallel per-binB cache of B-side payloads (one entry per fetched B
  # record, keyed by its source bOffset). Populated during fetchTile,
  # consumed when emitting each pair. Bounded by job size.
  var bExtras: Table[int, Table[int64, B]]
  for binB in sortedBinsB:
    bExtras[binB] = initTable[int64, B]()

  # Scratch buffers for INFO field decodes.
  var endData, svlenData, boffData: seq[int32]

  # Region-query one tile from the per-binB slim BCF. Lazily opens the reader
  # on first call for each binB. CSI queries return records whose [POS, END)
  # overlaps the region, so a record straddling a tile boundary appears in
  # both adjacent fetches; filter by POS to assign each record to exactly one
  # tile.
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
      let bOff = readBoff(vb, boffData)
      result.add(BufferedRec(
        pos: vb.POS, endPos: endB, id: $vb.ID, bOffset: bOff,
      ))
      bExtras[binB][bOff] = extract(vb)

  for va in vcfA.query(job.chrom):
    var endA: int64
    if not extractEnd(va, endData, svlenData, endA): continue
    let posA = va.POS
    let aOff = readBoff(va, boffData)

    for binB in sortedBinsB:
      let b = binB  # owned copy — lent iterator vars can't be captured
      let cands = buffers[b].getCandidates(posA, endA,
        proc(ti: int): seq[BufferedRec] = fetchTile(b, ti))
      for cand in cands:
        let sim =
          if cfg.metric == mOverlap:
            reciprocalOverlap(posA, endA, cand.pos, cand.endPos)
          else:
            jaccard(posA, endA, cand.pos, cand.endPos)
        if sim < cfg.threshold: continue
        let bExtra =
          if cand.bOffset in bExtras[b]: bExtras[b][cand.bOffset]
          else: default(B)
        let pr = emit(va, posA, endA, aOff, cand, bExtra, sim)
        if pr.keep:
          result.add(pr.item)

    # Evict tiles no future A record can need; payloads for those B records
    # become unreachable too — drop them in parallel to keep memory bounded.
    for binB, buf in buffers.mpairs:
      buf.evict(posA)

  vcfA.close()
  for v in vcfsB.mvalues:
    v.close()

# ---------------------------------------------------------------------------
# BND streaming (point events; slop-based proximity)
# ---------------------------------------------------------------------------

type
  BndPairEmitCb*[B, R] = proc(va: Variant; posA, pos2A, aOff: int64;
                              cand: BufferedRec; bExtra: B;
                              sim: float64): PairResult[R] {.closure.}

proc streamBndJobPairs*[B, R](job: MatchJob; cfg: MatchConfig;
                              extract: BExtractCb[B];
                              emit: BndPairEmitCb[B, R]): seq[R] =
  ## Per-job BND matching with a sliding cache of B records.
  ##
  ## A is streamed in POS order from the per-chrom sorted slim BCF, so the
  ## per-A window (posA - slop, posA + slop) advances monotonically. We
  ## maintain a Deque of B records currently in the active window, and on
  ## each A advance:
  ##   1. Evict from the front records whose POS dropped out of the window.
  ##   2. CSI-query *only* the new right-edge slice (cacheEnd .. winHi) we
  ##      haven't fetched yet, and append the decoded records to the cache.
  ##   3. Iterate the cache to filter by chr2/pos2 and dispatch `emit`.
  ##
  ## Net effect: each B record is fetched, decoded, and `extract()`-ed at
  ## most once across all A records whose windows include it. When A jumps
  ## past the cache (gap > 2*slop), the cache empties via eviction and the
  ## next delta query covers the full window — same cost as the old code.
  let slop = cfg.bndSlop
  if slop <= 0: return
  if 0 notin job.binsB: return            # B has no BND temp BCF

  var vcfA: VCF
  if not open(vcfA, job.pathA):
    raise newException(IOError, "cannot open A BND BCF: " & job.pathA)
  var vcfB: VCF
  if not open(vcfB, job.binsB[0]):
    raise newException(IOError, "cannot open B BND BCF: " & job.binsB[0])

  type BndCacheRec = object
    pos:    int64
    pos2:   int64
    bOff:   int64
    chr2:   string
    id:     string
    extra:  B

  var cache = initDeque[BndCacheRec]()
  ## Exclusive upper bound of the POS range currently covered by the cache
  ## (or by the union of all queries issued so far on this chrom). Records
  ## with pos in [cache.peekFirst().pos, cacheEnd) are guaranteed to be
  ## either in `cache` or to have been evicted because no future A can
  ## reach them. low(int64) forces a full first-A query.
  var cacheEnd: int64 = low(int64)

  var boffData, pos2Data: seq[int32]
  var chr2Data: string
  let twoSlop = float64(2 * slop)

  for va in vcfA.query(job.chrom):
    let posA = va.POS
    let aOff = readBoff(va, boffData)
    let p2A = readPos2(va, pos2Data)
    if not p2A.ok: continue
    let c2A = readChr2(va, chr2Data)
    if not c2A.ok: continue
    let pos2A = p2A.pos2
    let chr2A = c2A.chr2           # value copy; chr2Data scratch is reused below

    # Strict-band window (open both ends): pos must satisfy
    #   posA - slop < pos < posA + slop
    # i.e. inclusive integer range [posA - slop + 1, posA + slop - 1].
    let winLo = posA - slop.int64 + 1
    let winHi = posA + slop.int64           # exclusive upper bound

    # Evict from the left: records with pos < winLo can't match any current
    # or future A (A advances monotonically).
    while cache.len > 0 and cache.peekFirst().pos < winLo:
      cache.popFirst()

    # Delta query: fetch only the part of the window we don't already cover.
    let queryLo = max(cacheEnd, winLo)
    let queryHi = winHi
    if queryLo < queryHi:
      let regLo = max(1'i64, queryLo)
      let regHi = queryHi - 1                ## inclusive upper bound for hts region
      if regHi >= regLo:
        let region = job.chrom & ":" & $regLo & "-" & $regHi
        for vb in vcfB.query(region):
          let posB = vb.POS
          # Skip records the previous query already cached (boundary dedup).
          if posB < queryLo: continue
          let p2B = readPos2(vb, pos2Data)
          if not p2B.ok: continue
          let c2B = readChr2(vb, chr2Data)
          if not c2B.ok: continue
          cache.addLast(BndCacheRec(
            pos: posB, pos2: p2B.pos2, bOff: readBoff(vb, boffData),
            chr2: c2B.chr2, id: $vb.ID, extra: extract(vb),
          ))
      cacheEnd = queryHi

    # Emit: scan the cache (now == the window contents) and dispatch
    # passing pairs. Indexed iteration avoids a generic-deque `items`
    # resolution issue when BndCacheRec depends on the proc's type
    # parameter `B`.
    for i in 0 ..< cache.len:
      let cand = cache[i]
      if cand.chr2 != chr2A: continue
      let d2 = abs(cand.pos2 - pos2A)
      if d2 >= slop: continue
      let d1 = abs(cand.pos - posA)
      let sim = (twoSlop - float64(d1) - float64(d2)) / twoSlop
      if sim <= 0: continue                  ## safety net (band check enforces)
      let buffered = BufferedRec(
        pos: cand.pos, endPos: cand.pos + 1, id: cand.id, bOffset: cand.bOff,
      )
      let pr = emit(va, posA, pos2A, aOff, buffered, cand.extra, sim)
      if pr.keep:
        result.add(pr.item)

  vcfA.close()
  vcfB.close()
