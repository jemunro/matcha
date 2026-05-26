## Tests for src/matcha/bins.nim — bin math and tiled buffer.
## Run from project root: nim c --hints:off -r tests/test_bins.nim
echo "--------------- Test Bins ---------------"

import std/sets
import test_utils
import matcha/bins
import matcha/utils

# ---------------------------------------------------------------------------
# Bin assignment / range
# ---------------------------------------------------------------------------

# B01 — bin 0 catches everything < 512
timed("B01", "binIndexFor: SVLEN < 512 → bin 0"):
  doAssert binIndexFor(1) == 0
  doAssert binIndexFor(511) == 0

# B02 — bin 1 starts at 512
timed("B02", "binIndexFor: SVLEN 512 → bin 1, 1023 → bin 1"):
  doAssert binIndexFor(512) == 1
  doAssert binIndexFor(1023) == 1

# B03 — bin 2 starts at 1024
timed("B03", "binIndexFor: SVLEN 1024 → bin 2, 2047 → bin 2"):
  doAssert binIndexFor(1024) == 2
  doAssert binIndexFor(2047) == 2

# B04 — large SV
timed("B04", "binIndexFor: SVLEN 100000 → bin 8 ([2^16=65536, 2^17=131072))"):
  doAssert binIndexFor(100000) == 8

# B05 — degenerate values clamp to bin 0
timed("B05", "binIndexFor: 0 and negative clamp to bin 0"):
  doAssert binIndexFor(0) == 0
  doAssert binIndexFor(-1) == 0

# B06 — round-trip: binIndexFor(binRange(i).lo) == i, .hi-1 == i
timed("B06", "binRange round-trip for i=0,1,5,10"):
  for i in [0, 1, 5, 10]:
    let (lo, hi) = binRange(i)
    doAssert binIndexFor(lo) == i,
      "binIndexFor(binRange(" & $i & ").lo=" & $lo & ") != " & $i
    doAssert binIndexFor(hi - 1) == i,
      "binIndexFor(binRange(" & $i & ").hi-1=" & $(hi - 1) & ") != " & $i

# B07 — adjacent bins for binA=2 at threshold 0.5
timed("B07", "adjacentBins(2, 0.5, populated={0..5}) → {1,2,3}"):
  # bin 2 spans [1024, 2048); allowed B sizes (512, 4096) → bins 1,2,3.
  # Bin 0 [0,512): hi=512, qmin=1024*0.5=512. bHi > qmin? 512 > 512 false → out.
  # Bin 1 [512,1024): bLo=512 < qmax=4096, bHi=1024 > qmin=512 → in.
  # Bin 2 [1024,2048): in.
  # Bin 3 [2048,4096): bLo=2048 < 4096, bHi=4096 > 512 → in.
  # Bin 4 [4096,8192): bLo=4096 < 4096? false → out.
  let pop: set[uint8] = {0'u8, 1'u8, 2'u8, 3'u8, 4'u8, 5'u8}
  let adj = adjacentBins(2, 0.5, pop)
  let adjSet = adj.toHashSet
  doAssert adjSet == [1, 2, 3].toHashSet,
    "expected {1,2,3}, got " & $adj

# B08 — adjacent bins exclude populated bins outside the range
timed("B08", "adjacentBins skips populated bins outside [L_lo*t, L_hi/t)"):
  # binA=2 at t=0.5: window (512, 4096). Populated={0,5} → {} (both outside).
  let adj = adjacentBins(2, 0.5, {0'u8, 5'u8})
  doAssert adj.len == 0, "expected empty, got " & $adj

# B09 — t=0.9: only the same bin and immediate neighbours
timed("B09", "adjacentBins(2, 0.9): tight window picks bins 1,2,3"):
  # binA=2 [1024,2048); window (1024*0.9=921.6, 2048/0.9≈2275.6).
  # Bin 1 [512,1024): bLo=512<2275.6 ✓, bHi=1024>921.6 ✓ → in.
  # Bin 2 [1024,2048): in.
  # Bin 3 [2048,4096): bLo=2048<2275.6 ✓, bHi=4096>921.6 ✓ → in.
  let adj = adjacentBins(2, 0.9, {0'u8, 1'u8, 2'u8, 3'u8, 4'u8})
  let adjSet = adj.toHashSet
  doAssert adjSet == [1, 2, 3].toHashSet, "got " & $adj

# ---------------------------------------------------------------------------
# TiledBuffer
# ---------------------------------------------------------------------------

proc makeRec(pos, endPos: int64): BufferedRec =
  BufferedRec(pos: pos, endPos: endPos, srcIndex: 0, fileIdx: 0)

# T01 — cold buffer fetch loads expected tiles
timed("T01", "TiledBuffer: cold fetch loads tiles in expected range"):
  var buf = initTiledBuffer(1000, "chr1")
  var calls: seq[int]
  proc fetch(idx: int): seq[BufferedRec] =
    calls.add(idx)
    @[makeRec(idx.int64 * 1000 + 100, idx.int64 * 1000 + 500)]
  # posA=2500, queryEnd=2700. queryStart=1500. firstTile=1, lastTile=2.
  let cands = buf.getCandidates(2500, 2700, fetch)
  doAssert calls == @[1, 2], "expected fetch tiles [1,2], got " & $calls
  # Tile 1 has rec at (1100, 1500). pos=1100 < 2700, endPos=1500 > 2500? false → excluded.
  # Tile 2 has rec at (2100, 2500). endPos=2500 > 2500? false → excluded.
  doAssert cands.len == 0, "no records overlap [2500,2700)"

# T02 — re-query in same range: no extra fetch
timed("T02", "TiledBuffer: re-query reuses cached tiles"):
  var buf = initTiledBuffer(1000, "chr1")
  var calls: seq[int]
  proc fetch(idx: int): seq[BufferedRec] =
    calls.add(idx)
    @[makeRec(idx.int64 * 1000 + 200, idx.int64 * 1000 + 800)]
  discard buf.getCandidates(1500, 1700, fetch)  # tile 0,1
  let calls1 = calls.len
  discard buf.getCandidates(1500, 1700, fetch)  # same tiles
  doAssert calls.len == calls1, "expected no new fetch, got " & $(calls.len - calls1)

# T03 — empty tile recorded as fetched, no re-query
timed("T03", "TiledBuffer: empty fetch result still marks tile fetched"):
  var buf = initTiledBuffer(1000, "chr1")
  var fetchCount = 0
  proc fetch(idx: int): seq[BufferedRec] =
    inc fetchCount
    @[]
  discard buf.getCandidates(500, 700, fetch)  # tile 0 only (queryStart=0, queryEnd=700)
  doAssert fetchCount == 1
  discard buf.getCandidates(500, 700, fetch)
  doAssert fetchCount == 1, "should not re-fetch empty tile"

# T04 — A record spanning two tiles loads both
timed("T04", "TiledBuffer: span across tile boundary loads both"):
  var buf = initTiledBuffer(1000, "chr1")
  var fetched: seq[int]
  proc fetch(idx: int): seq[BufferedRec] =
    fetched.add(idx)
    @[]
  # posA=1900, queryEnd=2200. queryStart=900. firstTile=0, lastTile=2.
  discard buf.getCandidates(1900, 2200, fetch)
  doAssert fetched == @[0, 1, 2], "expected tiles [0,1,2], got " & $fetched

# T05 — evict drops only safely-old tiles
timed("T05", "TiledBuffer: evict preserves current and previous tile"):
  var buf = initTiledBuffer(1000, "chr1")
  proc fetch(idx: int): seq[BufferedRec] = @[]
  # Force tiles 0..5 into fetched
  for posA in [500'i64, 1500, 2500, 3500, 4500, 5500]:
    discard buf.getCandidates(posA, posA + 100, fetch)
  # posA=5500 → limit = floor(5500/1000) - 1 = 4. Evict tiles k<4 → {0,1,2,3}.
  buf.evict(5500)
  doAssert 4 in buf.fetched, "tile 4 should still be present"
  doAssert 5 in buf.fetched, "tile 5 should still be present"
  doAssert 3 notin buf.fetched, "tile 3 should be evicted"
  doAssert 0 notin buf.fetched, "tile 0 should be evicted"

# ---------------------------------------------------------------------------
# signedSvlen — VCF-spec sign convention applied at output write time
# ---------------------------------------------------------------------------

# S01 — DEL emits negative SVLEN
timed("S01", "signedSvlen: DEL → negative"):
  doAssert signedSvlen(svDEL, 100) == -100'i32
  doAssert signedSvlen(svDEL, 1)   ==   -1'i32

# S02 — DUP/INS/INV emit positive SVLEN
timed("S02", "signedSvlen: DUP/INS/INV → positive"):
  doAssert signedSvlen(svDUP, 100) == 100'i32
  doAssert signedSvlen(svINS, 100) == 100'i32
  doAssert signedSvlen(svINV, 100) == 100'i32

# S03 — BND and unknown types emit 0
timed("S03", "signedSvlen: BND/TRA/UNKNOWN → 0"):
  doAssert signedSvlen(svBND, 100)     == 0'i32
  doAssert signedSvlen(svTRA, 100)     == 0'i32
  doAssert signedSvlen(svUNKNOWN, 100) == 0'i32

# S04 — magnitude 0 stays 0 regardless of type
timed("S04", "signedSvlen: 0 magnitude → 0 for any svtype"):
  doAssert signedSvlen(svDEL, 0) == 0'i32
  doAssert signedSvlen(svDUP, 0) == 0'i32
  doAssert signedSvlen(svINS, 0) == 0'i32
