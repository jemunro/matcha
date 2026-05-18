## intervals.nim — pure interval overlap math for SV matching.
##
## Interval convention: [POS, END) half-open where length = END - POS.
## END is the value from INFO/END (1-based in VCF, treated here as exclusive
## upper bound so that length = END - POS = abs(SVLEN) for canonical SVs).

proc overlapData(posA, endA, posB, endB: int64): tuple[lenA, lenB, overlap: int64] =
  let lenA = endA - posA
  let lenB = endB - posB
  if lenA <= 0 or lenB <= 0: return
  let ov = min(endA, endB) - max(posA, posB)
  result = (lenA, lenB, max(0, ov))

proc reciprocalOverlap*(posA, endA, posB, endB: int64): float64 =
  ## Reciprocal overlap = overlap / max(lenA, lenB).
  ## Equivalent to requiring overlap/lenA >= t AND overlap/lenB >= t
  ## simultaneously (the standard definition used by truvari/bedtools).
  ## Returns 0.0 if intervals do not overlap or either has degenerate length.
  let (lenA, lenB, overlap) = overlapData(posA, endA, posB, endB)
  if overlap <= 0: return 0.0
  result = float64(overlap) / float64(max(lenA, lenB))

proc jaccard*(posA, endA, posB, endB: int64): float64 =
  ## Jaccard index = overlap / union.
  ## Returns 0.0 if intervals do not overlap or either has degenerate length.
  let (_, _, overlap) = overlapData(posA, endA, posB, endB)
  if overlap <= 0: return 0.0
  let unionLen = max(endA, endB) - min(posA, posB)
  if unionLen <= 0: return 0.0
  result = float64(overlap) / float64(unionLen)
