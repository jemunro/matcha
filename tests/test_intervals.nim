## Tests for src/matcha/intervals.nim — pure interval math.
## Run from project root: nim c --hints:off -r tests/test_intervals.nim
echo "--------------- Test Intervals ---------------"

import std/math
import test_utils
import matcha/intervals

# I01 — reciprocalOverlap: exact match
timed("I01", "reciprocalOverlap: exact match returns 1.0"):
  doAssert reciprocalOverlap(1000, 2000, 1000, 2000) == 1.0

# I02 — reciprocalOverlap: no overlap (gap)
timed("I02", "reciprocalOverlap: no overlap returns 0.0"):
  doAssert reciprocalOverlap(1000, 2000, 3000, 4000) == 0.0

# I03 — reciprocalOverlap: partial overlap correct value
timed("I03", "reciprocalOverlap: partial overlap A=[1000,2000) B=[1100,2100) = 0.9"):
  # overlap=900, min_len=1000 → 0.9
  let ro = reciprocalOverlap(1000, 2000, 1100, 2100)
  doAssert abs(ro - 0.9) < 1e-9, "expected 0.9, got " & $ro

# I04 — reciprocalOverlap: adjacent intervals (touching at boundary) = 0.0
timed("I04", "reciprocalOverlap: adjacent intervals return 0.0"):
  doAssert reciprocalOverlap(1000, 2000, 2000, 3000) == 0.0

# I05 — reciprocalOverlap: degenerate interval (len=0) = 0.0
timed("I05", "reciprocalOverlap: degenerate interval (posA==endA) returns 0.0"):
  doAssert reciprocalOverlap(1000, 1000, 1000, 2000) == 0.0

# I06 — jaccard: exact match
timed("I06", "jaccard: exact match returns 1.0"):
  doAssert jaccard(1000, 2000, 1000, 2000) == 1.0

# I07 — jaccard: no overlap
timed("I07", "jaccard: no overlap returns 0.0"):
  doAssert jaccard(1000, 2000, 3000, 4000) == 0.0

# I08 — jaccard: partial overlap
timed("I08", "jaccard: partial overlap A=[1000,2000) B=[1100,2100) = 900/1100"):
  # overlap=900, union=max(2000,2100)-min(1000,1100)=1100 → 900/1100
  let jac = jaccard(1000, 2000, 1100, 2100)
  doAssert abs(jac - 900.0 / 1100.0) < 1e-9, "expected " & $(900.0/1100.0) & ", got " & $jac

# I09 — jaccard: adjacent intervals = 0.0
timed("I09", "jaccard: adjacent intervals return 0.0"):
  doAssert jaccard(1000, 2000, 2000, 3000) == 0.0

# I10 — reciprocalOverlap: size asymmetry — large A, small B fully inside A
timed("I10", "reciprocalOverlap: large A=[17000,22000) small B=[17000,18000) = 0.2"):
  # overlap=1000, max_len=max(5000,1000)=5000 → 0.2
  let ro = reciprocalOverlap(17000, 22000, 17000, 18000)
  doAssert abs(ro - 0.2) < 1e-9, "expected 0.2, got " & $ro

# I11 — jaccard: size asymmetry — large A, small B
timed("I11", "jaccard: large A=[17000,22000) small B=[17000,18000) = 1000/5000 = 0.2"):
  # overlap=1000, union=5000 → 0.2
  let jac = jaccard(17000, 22000, 17000, 18000)
  doAssert abs(jac - 0.2) < 1e-9, "expected 0.2, got " & $jac

# I12 — queryWindow: threshold=0.5, svlen=1000
timed("I12", "queryWindow: svlen=1000, threshold=0.5 returns 500"):
  doAssert queryWindow(1000, 0.5) == 500

# I13 — queryWindow: threshold=0.0 returns svlen (full span)
timed("I13", "queryWindow: threshold=0.0 returns svlen"):
  doAssert queryWindow(1000, 0.0) == 1000

# I14 — queryWindow: threshold=1.0 returns 0
timed("I14", "queryWindow: threshold=1.0 returns 0"):
  doAssert queryWindow(1000, 1.0) == 0

# I15 — queryWindow: svlen=0 returns 0
timed("I15", "queryWindow: svlen=0 returns 0"):
  doAssert queryWindow(0, 0.5) == 0
