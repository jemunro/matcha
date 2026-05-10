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
    chrom*:    string
    posA*:     int64
    endA*:     int64
    idA*:      string
    posB*:     int64
    endB*:     int64
    idB*:      string
    svtype*:   SvType
    overlap*:  float64
    jaccard*:  float64
    aOffset*:  int64    ## BGZF virtual offset of A's source record (0 if unknown).
                        ## Populated for milestone-2 anno/merge/collapse;
                        ## not emitted in match TSV output.
    bOffset*:  int64    ## Same as aOffset, for B.

  MatchConfig* = object
    minOverlap*:      float64
    minJaccard*:      float64
    minOverlapSet*:   bool
    minJaccardSet*:   bool
    nThreads*:        int
    tmpDir*:          string
    outputPath*:      string
    callsetA*:        string
    callsetB*:        string
    selfMode*:        bool   ## When true, match callsetA against itself.
                             ## callsetB is empty; pair dedup uses MATCHA_BOFF.

const SupportedSvTypes* = {svDEL, svDUP, svINV}

const OutputHeader* =
  "#CHROM\tPOS_A\tEND_A\tID_A\tPOS_B\tEND_B\tID_B\tSVTYPE\tOVERLAP\tJACCARD"

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
  r.chrom & "\t" & $r.posA & "\t" & $r.endA & "\t" & r.idA & "\t" &
  $r.posB & "\t" & $r.endB & "\t" & r.idB & "\t" & $r.svtype & "\t" &
  formatFloat(r.overlap, ffDecimal, 6) & "\t" &
  formatFloat(r.jaccard, ffDecimal, 6)
