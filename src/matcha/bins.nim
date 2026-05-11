## bins.nim — log2 size-bin assignment and tiled B-record buffer.
##
## Size bins partition SV lengths on a log2 scale starting at 1024bp.
## Bin 0 catches everything < 1024bp; bin N (N >= 1) spans [2^(N+9), 2^(N+10)).
##
## Adjacent-bin pruning: for canonical reciprocal overlap (overlap/max) >= t,
## a B record of length L_b can only match an A record of length L_a if
##   L_b / L_a >= t  and  L_a / L_b >= t
## i.e. L_b ∈ [L_a * t, L_a / t].
## Given A is in bin [L_lo, L_hi), adjacent B bins are those whose range
## intersects (L_lo * t, L_hi / t).
##
## TiledBuffer: caches B records in fixed-width tiles (width = bin upper),
## loads lazily, and evicts tiles once A has advanced far enough that no
## future A record can need them.

import std/[bitops, sequtils, sets, tables]

const BinZeroUpper = 1024'i64

proc binIndexFor*(svlen: int64): int =
  ## Assign a bin index to a normalised (positive) SVLEN.
  ## Bin 0 = [0, 1024). Bin N = [2^(N+9), 2^(N+10)) for N >= 1.
  if svlen <= 0 or svlen < BinZeroUpper:
    return 0
  int(fastLog2(uint64(svlen))) - 9

proc binRange*(idx: int): tuple[lo, hi: int64] =
  ## Return the half-open [lo, hi) size range for a bin index.
  if idx <= 0:
    return (0'i64, BinZeroUpper)
  (1'i64 shl (idx + 9), 1'i64 shl (idx + 10))

proc adjacentBins*(binA: int, threshold: float64,
                   populatedB: set[uint8]): seq[int] =
  ## Return B bin indexes from populatedB whose size range intersects
  ## the set of B sizes that can satisfy reciprocalOverlap(A,B) >= threshold.
  ## For bin A with range [L_lo, L_hi), eligible B sizes are (L_lo*t, L_hi/t).
  ##
  ## Note: with overlap/max reciprocal overlap, both intervals must cover at
  ## least fraction t of the larger, bounding the size ratio to [t, 1/t].
  let (lo, hi) = binRange(binA)
  let qmin = float64(lo) * threshold
  let qmax = if threshold > 0.0: float64(hi) / threshold else: float64(high(int64))
  for i in populatedB:
    let (bLo, bHi) = binRange(int(i))
    if float64(bLo) < qmax and float64(bHi) > qmin:
      result.add(int(i))

# ---------------------------------------------------------------------------
# TiledBuffer — per-(svtype, binB) cache of B records
# ---------------------------------------------------------------------------

type
  BufferedRec* = object
    pos*:     int64
    endPos*:  int64
    id*:      string
    bOffset*: int64   ## BGZF virtual offset into the original source file
    chr2*:    string  ## mate chromosome (svBND only; empty for interval records)
    pos2*:    int64   ## mate position   (svBND only; zero  for interval records)

  TiledBuffer* = object
    tileWidth*: int64
    chrom*:     string
    tiles*:     Table[int, seq[BufferedRec]]
    fetched*:   HashSet[int]   ## tile indexes already attempted (even if empty)

proc initTiledBuffer*(tileWidth: int64, chrom: string): TiledBuffer =
  TiledBuffer(tileWidth: tileWidth, chrom: chrom)

proc getCandidates*(buf: var TiledBuffer, posA, queryEnd: int64,
                    fetchTile: proc(tileIdx: int): seq[BufferedRec]
                   ): seq[BufferedRec] =
  ## Return all BufferedRecs that could overlap [posA, queryEnd) from
  ## the B bin this buffer tracks. Loads tiles lazily via fetchTile.
  ##
  ## queryStart (left edge of the query window) = posA - tileWidth.
  ## tileWidth equals the upper bound of this bin's size range, which is
  ## the maximum SVLEN any cached B record can have — so a B record
  ## starting earlier than posA - tileWidth cannot extend far enough
  ## rightward to reach posA.
  if buf.tileWidth <= 0 or queryEnd <= posA:
    return

  let queryStart = max(0'i64, posA - buf.tileWidth)
  let firstTile  = int(queryStart div buf.tileWidth)
  let lastTile   = int((queryEnd - 1) div buf.tileWidth)

  for tileIdx in firstTile .. lastTile:
    if tileIdx notin buf.fetched:
      buf.tiles[tileIdx] = fetchTile(tileIdx)
      buf.fetched.incl(tileIdx)
    if tileIdx in buf.tiles:
      for r in buf.tiles[tileIdx]:
        if r.pos < queryEnd and r.endPos > posA:
          result.add(r)

proc evict*(buf: var TiledBuffer, posA: int64) =
  ## Discard tiles that no future A record (posA' >= posA) can need.
  ## A' needs tile K iff posA' < (K+2)*W  (because queryStart = posA' - W,
  ## and tile K covers up to (K+1)*W). Evict K when posA >= (K+2)*W,
  ## i.e. K < floor(posA/W) - 1.
  if buf.tileWidth <= 0:
    return
  let limit = int(posA div buf.tileWidth) - 1
  for k in toSeq(buf.fetched):
    if k < limit:
      buf.fetched.excl(k)
      buf.tiles.del(k)
