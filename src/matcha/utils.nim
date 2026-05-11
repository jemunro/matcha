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

  MatchResult* = object
    chrom*:      string
    posA*:       int64
    endA*:       int64    ## ignored for svBND rows; emitted as "." in output
    idA*:        string
    posB*:       int64
    endB*:       int64    ## ignored for svBND rows; emitted as "." in output
    idB*:        string
    svtype*:     SvType
    similarity*: float64
    aOffset*:    int64    ## BGZF virtual offset of A's source record (0 if unknown).
                          ## Populated for milestone-2 anno/merge/collapse;
                          ## not emitted in match TSV output.
    bOffset*:    int64    ## Same as aOffset, for B.

  MatchConfig* = object
    minOverlap*:      float64
    minJaccard*:      float64
    minOverlapSet*:   bool
    minJaccardSet*:   bool
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
  "#CHROM\tPOS_A\tEND_A\tID_A\tPOS_B\tEND_B\tID_B\tSVTYPE\tSIMILARITY"

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
  r.chrom & "\t" & $r.posA & "\t" & endAStr & "\t" & r.idA & "\t" &
  $r.posB & "\t" & endBStr & "\t" & r.idB & "\t" & $r.svtype & "\t" &
  formatFloat(r.similarity, ffDecimal, 6)
