## utils.nim — shared types and helpers for matcha.

import std/strutils

const NO_MATCH* = -1'i32  ## Sentinel for srcIndexB / fileIdxB in singleton MatchPairs.

type
  SvType* = enum
    svDEL = "DEL"
    svDUP = "DUP"
    svINV = "INV"
    svBND = "BND"
    svINS = "INS"      ## out of scope — silently skipped in preproc
    svTRA = "TRA"      ## out of scope — warned and skipped in preproc
    svUNKNOWN = "UNKNOWN"

  Metric* = enum
    ## Active interval-match metric. Exactly one is chosen per run via
    ## `--min-overlap` or `--min-jaccard`. BND matching is independent of
    ## this — it uses `--bnd-slop`.
    mOverlap = "overlap"
    mJaccard = "jaccard"

  MatchPair* = object
    ## 28-byte match record produced by matchcore.
    ## Layout: 4+4+4+4+4+2+2+2+1+1(pad) = 28 bytes.
    ## posA/posB enable CSI `chrom:pos-pos` queries for O(1) slim-BCF resolution;
    ## SRC_INDEX is the tiebreaker when multiple SVs share a position.
    srcIndexA*: int32   ## Sequential index of A record (identity + join key).
    srcIndexB*: int32   ## Sequential index of B record; NO_MATCH for singletons.
    posA*:      int32   ## POS of A record — for CSI resolution query.
    posB*:      int32   ## POS of B record; 0 for singletons.
    sim*:       float32 ## Similarity score in [0,1]; 0.0 for singletons.
    fileIdxA*:  int16   ## Index into the run's slim-BCF file list for A.
    fileIdxB*:  int16   ## Index into the run's slim-BCF file list for B; NO_MATCH for singletons.
    chromIdx*:  int16   ## Index into chromOrder (chrom shared by A and B in a job).
    svtype*:    int8    ## SvType cast to int8; read back as SvType(pair.svtype).

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
                             ## callsetB is empty; pair dedup uses srcIndexA < srcIndexB.
    emitSingletons*:  bool   ## When true, matchcore emits a MatchPair for every A
                             ## record with no passing B match (collapse only).

const SupportedSvTypes* = {svDEL, svDUP, svINV, svBND}

const OutputHeader* =
  "#CHROM_A\tPOS_A\tEND_A\tID_A\tCHROM_B\tPOS_B\tEND_B\tID_B\tSVTYPE\tSIMILARITY"

proc isStdoutPath*(p: string): bool =
  ## True when `p` refers to standard output (empty, "-", or "/dev/stdout").
  p == "" or p == "-" or p == "/dev/stdout"

proc parseSvType*(s: string): SvType =
  case s.toUpperAscii
  of "DEL": svDEL
  of "DUP": svDUP
  of "INV": svINV
  of "BND": svBND
  of "INS": svINS
  of "TRA": svTRA
  else: svUNKNOWN
