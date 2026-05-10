## intervals.nim — pure interval overlap and window computation for SV matching.
##
## Interval convention: [POS, END) half-open where length = END - POS.
## END is the value from INFO/END (1-based in VCF, treated here as exclusive
## upper bound so that length = END - POS = abs(SVLEN) for canonical SVs).

import std/math

proc reciprocalOverlap*(posA, endA, posB, endB: int64): float64 =
  ## Reciprocal overlap = overlap / max(lenA, lenB).
  ## Equivalent to requiring overlap/lenA >= t AND overlap/lenB >= t
  ## simultaneously (the standard definition used by truvari/bedtools).
  ## Returns 0.0 if intervals do not overlap or either has degenerate length.
  let lenA = endA - posA
  let lenB = endB - posB
  if lenA <= 0 or lenB <= 0:
    return 0.0
  let overlapStart = max(posA, posB)
  let overlapEnd   = min(endA, endB)
  let overlap = overlapEnd - overlapStart
  if overlap <= 0:
    return 0.0
  result = float64(overlap) / float64(max(lenA, lenB))

proc jaccard*(posA, endA, posB, endB: int64): float64 =
  ## Jaccard index = overlap / union.
  ## Returns 0.0 if intervals do not overlap or either has degenerate length.
  let lenA = endA - posA
  let lenB = endB - posB
  if lenA <= 0 or lenB <= 0:
    return 0.0
  let overlapStart = max(posA, posB)
  let overlapEnd   = min(endA, endB)
  let overlap = overlapEnd - overlapStart
  if overlap <= 0:
    return 0.0
  let unionLen = max(endA, endB) - min(posA, posB)
  if unionLen <= 0:
    return 0.0
  result = float64(overlap) / float64(unionLen)

proc queryWindow*(svlen: int64, threshold: float64): int64 =
  ## Half-width of the candidate query window.
  ## window = ceil(svlen * (1.0 - threshold))
  ## svlen must be non-negative (callers pass abs(SVLEN)).
  ## Returns 0 if svlen <= 0 or threshold >= 1.0.
  if svlen <= 0 or threshold >= 1.0:
    return 0
  result = int64(ceil(float64(svlen) * (1.0 - threshold)))
