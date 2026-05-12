## utils.nim — shared types and helpers for matcha.

import std/strutils

type
  SvType* = enum
    svDEL = "DEL"
    svDUP = "DUP"
    svINV = "INV"
    svBND = "BND"
    svINS = "INS"
    svTRA = "TRA"
    svUNKNOWN = "UNKNOWN"

  Metric* = enum
    ## Active interval-match metric. Exactly one is chosen per run via
    ## `--min-overlap` or `--min-jaccard`. BND matching is independent of
    ## this — it uses `--bnd-slop`.
    mOverlap = "overlap"
    mJaccard = "jaccard"

  MatchResult* = object
    chromA*:     string
    posA*:       int64
    endA*:       int64    ## ignored for svBND rows; emitted as "." in output
    idA*:        string
    chromB*:     string
    posB*:       int64
    endB*:       int64    ## ignored for svBND rows; emitted as "." in output
    idB*:        string
    svtype*:     SvType
    similarity*: float64
    aOffset*:    int64    ## BGZF virtual offset of A's source record (0 if unknown).
                          ## Used by anno phase 3 to rejoin annotations against
                          ## the original input. Not emitted in TSV output.
    bOffset*:    int64    ## Same as aOffset, for B.

  MatchConfig* = object
    metric*:          Metric   ## Active interval metric (mOverlap | mJaccard).
    threshold*:       float64  ## Minimum score for the active metric (0.0-1.0).
    bndSlop*:         int      ## --bnd-slop (default 100); both breakends of a
                               ## BND pair must lie within this many bp.
    nThreads*:        int
    tmpDir*:          string
    outputPath*:      string
    callsetA*:        string
    callsetB*:        string
    selfMode*:        bool   ## When true, match callsetA against itself.
                             ## callsetB is empty; pair dedup uses MATCHA_BOFF.

const SupportedSvTypes* = {svDEL, svDUP, svINV, svBND}

const OutputHeader* =
  "#CHROM_A\tPOS_A\tEND_A\tID_A\tCHROM_B\tPOS_B\tEND_B\tID_B\tSVTYPE\tSIMILARITY"

proc parseSvType*(s: string): SvType =
  case s.toUpperAscii
  of "DEL": svDEL
  of "DUP": svDUP
  of "INV": svINV
  of "BND": svBND
  of "INS": svINS
  of "TRA": svTRA
  else: svUNKNOWN

proc formatMatchResult*(r: MatchResult): string =
  ## Format a MatchResult as a tab-separated line (no trailing newline).
  ## BND rows emit "." for END_A / END_B (BNDs are points, not intervals).
  let endAStr = if r.svtype == svBND: "." else: $r.endA
  let endBStr = if r.svtype == svBND: "." else: $r.endB
  r.chromA & "\t" & $r.posA & "\t" & endAStr & "\t" & r.idA & "\t" &
  r.chromB & "\t" & $r.posB & "\t" & endBStr & "\t" & r.idB & "\t" & $r.svtype & "\t" &
  formatFloat(r.similarity, ffDecimal, 6)
