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

import std/[algorithm, sequtils, tables]
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
                            ovl, jac: float64): PairResult[R] {.closure.}

proc streamJobPairs*[B, R](job: MatchJob; cfg: MatchConfig;
                           extract: BExtractCb[B];
                           emit: PairEmitCb[B, R]): seq[R] =
  ## Drive the per-job match: stream A from the per-(svtype, binA) BCF
  ## restricted to job.chrom, fetch B candidates lazily through TiledBuffers,
  ## compute overlap/jaccard, apply the threshold filter, and dispatch each
  ## passing pair to `emit`. The optional `extract` callback runs once per
  ## fetched B record and its return value is threaded into the matching
  ## `emit` call — used by anno mode to carry user-requested INFO values.
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
        let ovl = reciprocalOverlap(posA, endA, cand.pos, cand.endPos)
        let jac = jaccard(posA, endA, cand.pos, cand.endPos)
        let passOverlap = (not cfg.minOverlapSet) or (ovl >= cfg.minOverlap)
        let passJaccard = (not cfg.minJaccardSet) or (jac >= cfg.minJaccard)
        if passOverlap and passJaccard:
          let bExtra =
            if cand.bOffset in bExtras[b]: bExtras[b][cand.bOffset]
            else: default(B)
          let pr = emit(va, posA, endA, aOff, cand, bExtra, ovl, jac)
          if pr.keep:
            result.add(pr.item)

    # Evict tiles no future A record can need; payloads for those B records
    # become unreachable too — drop them in parallel to keep memory bounded.
    for binB, buf in buffers.mpairs:
      buf.evict(posA)

  vcfA.close()
  for v in vcfsB.mvalues:
    v.close()
