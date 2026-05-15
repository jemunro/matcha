## matchcore.nim — shared per-job matching loop.
##
## Both `matcha match`, `matcha anno`, and `matcha collapse` walk the same
## A-vs-B candidate space. This module is intentionally minimal: it returns
## pure (aOff, bOff, sim) triples and leaves per-mode field resolution
## (CHROM/POS/END/ID/INFO) to the adapters in match.nim and anno.nim. Those
## adapters re-scan the same slim BCFs that matchcore reads, picking records
## up by INFO/MATCHA_BOFF — the slim BCFs are tiny and already CSI-indexed.
##
## Self-mode dedup (aOff < bOff filter, dropping the trivial X-vs-X case) is
## baked in via cfg.selfMode; anno never sets it.

import std/[algorithm, deques, sequtils, sets, tables]
import hts
import intervals, preproc, bins, utils

type
  MatchPair* = object
    aOff*: int64
    bOff*: int64
    sim*:  float64

# ---------------------------------------------------------------------------
# Slim-BCF INFO decode helpers (also used by per-mode resolution adapters)
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
# Shared slim-BCF scan template
# ---------------------------------------------------------------------------

template scanSlimBcf*(path, chrom: string; needed: HashSet[int64];
                      body: untyped) =
  ## Open a slim BCF, query `chrom`, iterate records whose MATCHA_BOFF is in
  ## `needed`, and execute `body` for each matched record.
  ## Injects into `body`: `v` (Variant) and `off` (int64, the MATCHA_BOFF).
  block:
    var vcf: VCF
    if not open(vcf, path):
      raise newException(IOError, "cannot reopen slim BCF: " & path)
    var boffScratch: seq[int32]
    for v {.inject.} in vcf.query(chrom):
      let off {.inject.} = readBoff(v, boffScratch)
      if off notin needed: continue
      body
    vcf.close()

# ---------------------------------------------------------------------------
# Interval matching
# ---------------------------------------------------------------------------

proc streamJobPairs*(job: MatchJob; cfg: MatchConfig): seq[MatchPair] =
  ## Drive the per-job match: stream A from the per-(svtype, binA) BCF
  ## restricted to job.chrom, fetch B candidates lazily through TiledBuffers,
  ## compute the active interval metric, apply the threshold, and append a
  ## (aOff, bOff, sim) triple for each passing pair.
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

  var endData, svlenData, boffData: seq[int32]

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
      result.add(BufferedRec(
        pos: vb.POS, endPos: endB, bOffset: readBoff(vb, boffData),
      ))

  for va in vcfA.query(job.chrom):
    var endA: int64
    if not extractEnd(va, endData, svlenData, endA): continue
    let posA = va.POS
    let aOff = readBoff(va, boffData)

    for binB in sortedBinsB:
      let b = binB   # owned copy — lent iterator vars can't be captured
      let cands = buffers[b].getCandidates(posA, endA,
        proc(ti: int): seq[BufferedRec] = fetchTile(b, ti))
      for cand in cands:
        if cfg.selfMode and aOff >= cand.bOffset: continue
        let sim =
          if cfg.metric == mOverlap:
            reciprocalOverlap(posA, endA, cand.pos, cand.endPos)
          else:
            jaccard(posA, endA, cand.pos, cand.endPos)
        if sim < cfg.threshold: continue
        result.add(MatchPair(aOff: aOff, bOff: cand.bOffset, sim: sim))

    # Evict tiles no future A record can need.
    for binB, buf in buffers.mpairs:
      buf.evict(posA)

  vcfA.close()
  for v in vcfsB.mvalues:
    v.close()

# ---------------------------------------------------------------------------
# BND matching (point events; slop-based proximity)
# ---------------------------------------------------------------------------

proc streamBndJobPairs*(job: MatchJob; cfg: MatchConfig): seq[MatchPair] =
  ## Per-job BND matching with a sliding cache of B records.
  ##
  ## A is streamed in POS order from the per-chrom sorted slim BCF, so the
  ## per-A window (posA - slop, posA + slop) advances monotonically. We
  ## maintain a Deque of B records currently in the active window, and on
  ## each A advance:
  ##   1. Evict from the front records whose POS dropped out of the window.
  ##   2. CSI-query *only* the new right-edge slice (cacheEnd .. winHi) we
  ##      haven't fetched yet, and append the decoded records to the cache.
  ##   3. Iterate the cache to filter by chr2/pos2 and emit passing pairs.
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

  var cache = initDeque[BndCacheRec]()
  ## Exclusive upper bound of the POS range currently covered by the cache
  ## (or by the union of all queries issued so far on this chrom).
  ## low(int64) forces a full first-A query.
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
          if posB < queryLo: continue
          let p2B = readPos2(vb, pos2Data)
          if not p2B.ok: continue
          let c2B = readChr2(vb, chr2Data)
          if not c2B.ok: continue
          cache.addLast(BndCacheRec(
            pos: posB, pos2: p2B.pos2, bOff: readBoff(vb, boffData),
            chr2: c2B.chr2,
          ))
      cacheEnd = queryHi

    for i in 0 ..< cache.len:
      let cand = cache[i]
      if cand.chr2 != chr2A: continue
      let d2 = abs(cand.pos2 - pos2A)
      if d2 >= slop: continue
      let d1 = abs(cand.pos - posA)
      let sim = (twoSlop - float64(d1) - float64(d2)) / twoSlop
      if sim <= 0: continue                  ## safety net (band check enforces)
      if cfg.selfMode and aOff >= cand.bOff: continue
      result.add(MatchPair(aOff: aOff, bOff: cand.bOff, sim: sim))

  vcfA.close()
  vcfB.close()
