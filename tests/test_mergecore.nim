## Tests for src/matcha/mergecore.nim — clustering correctness.
## Run from project root: nim c --hints:off -r tests/test_mergecore.nim
echo "--------------- Test Mergecore ---------------"

import std/[algorithm, random, tables]
import test_utils
import matcha/mergecore

# Canonicalise a clustering for equality comparison: sort each cluster's
# members, then sort the cluster list by first member.
proc canonicalise(clusters: seq[seq[int32]]): seq[seq[int32]] =
  var out2 = newSeq[seq[int32]](clusters.len)
  for i, cl in clusters.pairs:
    var c = cl
    c.sort(system.cmp)
    out2[i] = c
  out2.sort(proc(a, b: seq[int32]): int = system.cmp(a[0], b[0]))
  out2

# Build a random simMap of N offsets [0..N-1] with given edge probability and
# similarity drawn uniformly from [0.0, 1.0]. Seeded for determinism.
proc randomSimMap(n: int; pEdge: float; seed: int64):
                  tuple[offsets: seq[int32], simMap: Table[(int32, int32), float64]] =
  var rng = initRand(seed)
  result.offsets = newSeq[int32](n)
  for i in 0 ..< n: result.offsets[i] = int32(i)
  for i in 0 ..< n:
    for j in i + 1 ..< n:
      if rng.rand(1.0) < pEdge:
        let sim = rng.rand(1.0)
        result.simMap[(int32(i), int32(j))] = sim

# C01 — agglomerateSparse matches agglomerateDense across linkages and thresholds
# at small N (full dense path).
timed("AC01", "dense vs sparse match: N=8, all linkages, 4 thresholds"):
  let (offsets, simMap) = randomSimMap(8, 0.6, 11)
  for lk in [lmSingle, lmComplete, lmAverage]:
    for th in [0.1, 0.3, 0.5, 0.8]:
      let d = canonicalise(agglomerateDense(offsets, simMap, lk, th))
      let s = canonicalise(agglomerateSparse(offsets, simMap, lk, th))
      doAssert d == s,
        "N=8 linkage=" & $lk & " threshold=" & $th & " differ:\n  dense=" & $d & "\n  sparse=" & $s

# AC02 — sparse-only-region N (above dense dispatcher threshold).
timed("AC02", "dense vs sparse match: N=50, dense edges, all linkages"):
  let (offsets, simMap) = randomSimMap(50, 0.4, 42)
  for lk in [lmSingle, lmComplete, lmAverage]:
    for th in [0.2, 0.5, 0.75]:
      let d = canonicalise(agglomerateDense(offsets, simMap, lk, th))
      let s = canonicalise(agglomerateSparse(offsets, simMap, lk, th))
      doAssert d == s,
        "N=50 linkage=" & $lk & " threshold=" & $th & " differ"

# AC03 — sparse graph: most pairs have no edge.
timed("AC03", "dense vs sparse match: N=200, sparse edges (p=0.05)"):
  let (offsets, simMap) = randomSimMap(200, 0.05, 7)
  for lk in [lmSingle, lmComplete, lmAverage]:
    for th in [0.3, 0.6]:
      let d = canonicalise(agglomerateDense(offsets, simMap, lk, th))
      let s = canonicalise(agglomerateSparse(offsets, simMap, lk, th))
      doAssert d == s,
        "N=200 sparse linkage=" & $lk & " threshold=" & $th & " differ"

# AC04 — beyond AggDenseThreshold: pure sparse path through the dispatcher,
# but we still compare against the dense reference for correctness.
timed("AC04", "dense vs sparse match: N=300, p=0.2"):
  let (offsets, simMap) = randomSimMap(300, 0.2, 99)
  for lk in [lmSingle, lmComplete, lmAverage]:
    let th = 0.5
    let d = canonicalise(agglomerateDense(offsets, simMap, lk, th))
    let s = canonicalise(agglomerateSparse(offsets, simMap, lk, th))
    doAssert d == s,
      "N=300 linkage=" & $lk & " threshold=" & $th & " differ"

# AC05 — degenerate: empty + singleton.
timed("AC05", "dense vs sparse: empty and singleton inputs"):
  var emptyMap: Table[(int32, int32), float64]
  doAssert agglomerateDense(@[], emptyMap, lmAverage, 0.5) == @[]
  doAssert agglomerateSparse(@[], emptyMap, lmAverage, 0.5) == @[]
  doAssert agglomerateDense(@[5'i32], emptyMap, lmAverage, 0.5) == @[@[5'i32]]
  doAssert agglomerateSparse(@[5'i32], emptyMap, lmAverage, 0.5) == @[@[5'i32]]

# AC06 — unreachable threshold: no edges can pass → all singletons.
timed("AC06", "no edges above threshold → singleton clusters"):
  let (offsets, simMap) = randomSimMap(20, 0.5, 3)
  for lk in [lmSingle, lmComplete, lmAverage]:
    let d = canonicalise(agglomerateDense(offsets, simMap, lk, 1.01))
    let s = canonicalise(agglomerateSparse(offsets, simMap, lk, 1.01))
    doAssert d == s
    doAssert d.len == 20  # 20 singletons

# AC07 — fully-connected component with identical similarities (tie-break
# stress test). Output must be deterministic and match between paths.
timed("AC07", "ties: all edges = same similarity"):
  var offsets: seq[int32]
  for i in 0 ..< 8: offsets.add(int32(i))
  var simMap: Table[(int32, int32), float64]
  for i in 0 ..< 8:
    for j in i + 1 ..< 8:
      simMap[(int32(i), int32(j))] = 0.7
  for lk in [lmSingle, lmComplete, lmAverage]:
    let d = canonicalise(agglomerateDense(offsets, simMap, lk, 0.5))
    let s = canonicalise(agglomerateSparse(offsets, simMap, lk, 0.5))
    doAssert d == s,
      "tie test linkage=" & $lk & " differ:\n  dense=" & $d & "\n  sparse=" & $s
