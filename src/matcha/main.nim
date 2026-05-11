## main.nim — CLI argument parsing and subcommand dispatch.
## Entry point is src/matcha.nim which includes this file.

import std/[os, parseopt, strutils]
import utils, match, anno, log

const NimblePkgVersion {.strdefine.} = "dev"
const VERSION = NimblePkgVersion

const ShortNoVal = {'h', 'v'}

proc nextVal(p: var OptParser; flag: string): string =
  ## Return the value for a flag, consuming the next argv token if needed.
  if p.val != "":
    return p.val
  p.next()
  if p.kind == cmdArgument:
    return p.key
  stderr.writeLine "error: --" & flag & " requires a value"
  quit(1)

proc usage(code: int = 1) =
  let f = if code == 0: stdout else: stderr
  f.writeLine "matcha v" & VERSION
  f.writeLine ""
  f.writeLine "Usage: matcha <subcommand> [options]"
  f.writeLine ""
  f.writeLine "Subcommands:"
  f.writeLine "  match    Find pairwise matches between two SV callsets"
  f.writeLine "  anno     Annotate a callset with INFO fields from a database VCF"
  f.writeLine ""
  f.writeLine "Flags:"
  f.writeLine "  --version       show version"
  f.writeLine "  -v, --verbose   verbose logging to stderr"
  f.writeLine "  -h, --help      show this help"
  f.writeLine ""
  f.writeLine "Run 'matcha <subcommand> --help' for subcommand options."
  quit(code)

proc matchUsage(code: int = 1) =
  let f = if code == 0: stdout else: stderr
  f.writeLine "Usage:"
  f.writeLine "  matcha match [options] callsetA callsetB   # match A against B"
  f.writeLine "  matcha match --self [options] INPUT        # match INPUT against itself"
  f.writeLine ""
  f.writeLine "Inputs may be VCF.gz (.vcf.gz) or BCF (.bcf); format is detected automatically."
  f.writeLine ""
  f.writeLine "Options:"
  f.writeLine "  --min-overlap FLOAT             minimum reciprocal overlap (0.0-1.0)"
  f.writeLine "  --min-jaccard FLOAT             minimum Jaccard index (0.0-1.0)"
  f.writeLine "  --bnd-slop INT                  max breakend offset for BND matches (default: 100)"
  f.writeLine "  --self                          match a single input against itself"
  f.writeLine "                                  (each pair emitted once; no self-self)"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  --output PATH                   output file (default: stdout)"
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      show this help"
  f.writeLine ""
  f.writeLine "Exactly one of --min-overlap or --min-jaccard is required."
  f.writeLine "BND rows are matched by --bnd-slop only; the chosen metric drives interval rows."
  f.writeLine ""
  f.writeLine "Output is tab-separated. A ##matcha_metric=<overlap|jaccard> preamble"
  f.writeLine "precedes the #-prefixed header line. Columns:"
  f.writeLine "  #CHROM  POS_A  END_A  ID_A  POS_B  END_B  ID_B  SVTYPE  SIMILARITY"
  f.writeLine "BND rows emit '.' for END_A / END_B."
  quit(code)

proc runMatch(rawArgs: seq[string]) =
  var cfg = MatchConfig(nThreads: 1, bndSlop: 100)
  var positionals: seq[string]
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      case p.key
      of "min-overlap":
        let v = nextVal(p, "min-overlap")
        try: cfg.minOverlap = parseFloat(v)
        except ValueError:
          stderr.writeLine "error: --min-overlap must be a float, got: " & v
          quit(1)
        cfg.minOverlapSet = true
      of "min-jaccard":
        let v = nextVal(p, "min-jaccard")
        try: cfg.minJaccard = parseFloat(v)
        except ValueError:
          stderr.writeLine "error: --min-jaccard must be a float, got: " & v
          quit(1)
        cfg.minJaccardSet = true
      of "bnd-slop":
        let v = nextVal(p, "bnd-slop")
        try: cfg.bndSlop = parseInt(v)
        except ValueError:
          stderr.writeLine "error: --bnd-slop must be an integer, got: " & v
          quit(1)
        if cfg.bndSlop <= 0:
          stderr.writeLine "error: --bnd-slop must be > 0"
          quit(1)
      of "threads":
        let v = nextVal(p, "threads")
        try: cfg.nThreads = parseInt(v)
        except ValueError:
          stderr.writeLine "error: --threads must be an integer, got: " & v
          quit(1)
        if cfg.nThreads < 1:
          stderr.writeLine "error: --threads must be >= 1"
          quit(1)
      of "tmp-dir":
        cfg.tmpDir = nextVal(p, "tmp-dir")
      of "output":
        cfg.outputPath = nextVal(p, "output")
      of "self":
        cfg.selfMode = true
      of "v", "verbose":
        setVerbose(true)
      of "h", "help":
        matchUsage(0)
      else:
        stderr.writeLine "error: unknown option: --" & p.key
        matchUsage()
    of cmdArgument:
      positionals.add(p.key)

  if cfg.minOverlapSet == cfg.minJaccardSet:
    if cfg.minOverlapSet:
      stderr.writeLine "error: --min-overlap and --min-jaccard are mutually exclusive"
    else:
      stderr.writeLine "error: exactly one of --min-overlap or --min-jaccard is required"
    matchUsage()
  let expected = if cfg.selfMode: 1 else: 2
  if positionals.len != expected:
    let what = if cfg.selfMode: "1 input file (--self mode)" else: "2 input files"
    stderr.writeLine "error: expected " & what & ", got " & $positionals.len
    matchUsage()

  cfg.callsetA = positionals[0]
  if not cfg.selfMode:
    cfg.callsetB = positionals[1]

  if not fileExists(cfg.callsetA):
    stderr.writeLine "error: input file not found: " & cfg.callsetA
    quit(1)
  if not cfg.selfMode and not fileExists(cfg.callsetB):
    stderr.writeLine "error: input file not found: " & cfg.callsetB
    quit(1)

  if cfg.tmpDir == "":
    cfg.tmpDir = getTempDir()

  runMatch(cfg)

proc annoUsage(code: int = 1) =
  ## Short usage shown on error paths. Lists synopsis, options, and one
  ## example call so the user can see how -a expressions are constructed
  ## without scrolling through the full docs. `--help` calls annoHelp() for
  ## the long-form aggregation-function and implicit-variable reference.
  let f = if code == 0: stdout else: stderr
  f.writeLine "Usage:"
  f.writeLine "  matcha anno [options] input database"
  f.writeLine ""
  f.writeLine "Options:"
  f.writeLine "  -a OUTFIELD=FUNC(SRCFIELD)      annotation expression (repeatable, >=1 required)"
  f.writeLine "  -o PATH                         output (.vcf | .vcf.gz | .bcf); default stdout VCF"
  f.writeLine "  --min-overlap FLOAT             minimum reciprocal overlap (0.0-1.0)"
  f.writeLine "  --min-jaccard FLOAT             minimum Jaccard index (0.0-1.0)"
  f.writeLine "  --bnd-slop INT                  max breakend offset for BND matches (default: 100)"
  f.writeLine "  --overwrite                     replace OUTFIELDs already in input header"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      show full help (functions, MATCHA_* variables)"
  f.writeLine ""
  f.writeLine "Exactly one of --min-overlap or --min-jaccard is required."
  f.writeLine ""
  f.writeLine "Example:"
  f.writeLine "  matcha anno --min-overlap 0.7 \\"
  f.writeLine "    -a AF_DB=max(AF) \\"
  f.writeLine "    -a CALLERS=unique(CALLERS) \\"
  f.writeLine "    -a N=first(MATCHA_COUNT) \\"
  f.writeLine "    -o annotated.vcf.gz input.vcf.gz gnomad_sv.vcf.gz"
  f.writeLine ""
  f.writeLine "Run 'matcha anno --help' for aggregation functions and MATCHA_* variables."
  quit(code)

proc annoHelp() =
  ## Full --help: synopsis + options + complete reference for annotation
  ## expressions, aggregation functions, and implicit MATCHA_* variables.
  let f = stdout
  f.writeLine "Usage:"
  f.writeLine "  matcha anno [options] input database"
  f.writeLine ""
  f.writeLine "Annotate input VCF/BCF with INFO fields aggregated from a database VCF/BCF,"
  f.writeLine "based on SV matches under the same coordinate+size criteria as `matcha match`."
  f.writeLine "Inputs may be VCF.gz (.vcf.gz) or BCF (.bcf); format is detected automatically."
  f.writeLine ""
  f.writeLine "Options:"
  f.writeLine "  -a OUTFIELD=FUNC(SRCFIELD)      annotation expression (repeatable, >=1 required)"
  f.writeLine "  -o PATH                         output file; format from extension:"
  f.writeLine "                                    .vcf     → uncompressed VCF"
  f.writeLine "                                    .vcf.gz  → bgzipped VCF (+ .csi index)"
  f.writeLine "                                    .bcf     → BCF (+ .csi index)"
  f.writeLine "                                  Default: uncompressed VCF to stdout."
  f.writeLine "  --min-overlap FLOAT             minimum reciprocal overlap (0.0-1.0)"
  f.writeLine "  --min-jaccard FLOAT             minimum Jaccard index (0.0-1.0)"
  f.writeLine "  --bnd-slop INT                  max breakend offset for BND matches (default: 100)"
  f.writeLine "  --overwrite                     replace OUTFIELDs that already exist in input header"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      show this help"
  f.writeLine ""
  f.writeLine "Exactly one of --min-overlap or --min-jaccard is required."
  f.writeLine ""
  f.writeLine "Annotation expressions"
  f.writeLine "  Form:  -a OUTFIELD=FUNC(SRCFIELD)"
  f.writeLine "    OUTFIELD   name of the INFO field to emit on the output."
  f.writeLine "               Must be unique across all -a expressions in a single run."
  f.writeLine "               If already present in the input header, --overwrite is required."
  f.writeLine "    FUNC       aggregation function applied across all matching database records"
  f.writeLine "               for each input record (see below)."
  f.writeLine "    SRCFIELD   an INFO field name in the database VCF, or an implicit MATCHA_*"
  f.writeLine "               variable (see below)."
  f.writeLine ""
  f.writeLine "  Examples:"
  f.writeLine "    -a AF_DB=max(AF)                  highest AF seen across matches"
  f.writeLine "    -a CALLERS=unique(CALLERS)        deduped flattened list of caller IDs"
  f.writeLine "    -a TOP_CALLER=best(CALLER)        CALLER on the best-scoring match"
  f.writeLine "    -a N_MATCH=first(MATCHA_COUNT)    number of database matches (0 if none)"
  f.writeLine "    -a TOP_SIM=max(MATCHA_SIMILARITY) highest similarity across matches"
  f.writeLine ""
  f.writeLine "Aggregation functions"
  f.writeLine "  max | min | mean    numeric.  Pool all values across matches, then aggregate."
  f.writeLine "                      'mean' always emits Float (integer source coerces)."
  f.writeLine "  first | last        any.       Value from the earliest / latest match by"
  f.writeLine "                                 database-record position."
  f.writeLine "  best                any.       Value from the match with the highest"
  f.writeLine "                                 SIMILARITY; ties broken by earliest position."
  f.writeLine "  all                 any.       All values, comma-separated. Output Number=."
  f.writeLine "  unique              any.       Deduplicated values, comma-separated."
  f.writeLine "                                 Output Number=."
  f.writeLine ""
  f.writeLine "  List-valued database INFO fields (Number=. or N>1) are flattened across all"
  f.writeLine "  matching records before the aggregation function runs."
  f.writeLine ""
  f.writeLine "Implicit MATCHA_* variables"
  f.writeLine "  Available as SRCFIELD in any -a expression without needing a database header"
  f.writeLine "  declaration:"
  f.writeLine "    MATCHA_COUNT      Integer (scalar).  Number of database matches for the"
  f.writeLine "                      input record. On unmatched records, expressions wrapping"
  f.writeLine "                      MATCHA_COUNT emit 0; other expressions leave their"
  f.writeLine "                      OUTFIELD absent."
  f.writeLine "    MATCHA_SIMILARITY Float (vector).  Per-match similarity. For interval"
  f.writeLine "                      matches this is the active metric (--min-overlap or"
  f.writeLine "                      --min-jaccard). For BND matches this is the slop-based"
  f.writeLine "                      proximity score (2*slop - |dPOS| - |dPOS2|) / (2*slop)."
  f.writeLine ""
  f.writeLine "  Note: best(MATCHA_SIMILARITY) and max(MATCHA_SIMILARITY) coincide, since"
  f.writeLine "  best() now ranks by similarity directly."
  quit(0)

proc runAnnoCli(rawArgs: seq[string]) =
  var cfg = AnnoConfig(nThreads: 1, bndSlop: 100)
  var positionals: seq[string]
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      case p.key
      of "a":
        let raw = nextVal(p, "a")
        try:
          cfg.exprs.add(parseAnnoExpr(raw))
        except ValueError as e:
          stderr.writeLine "error: " & e.msg
          quit(1)
      of "o", "output":
        cfg.outputPath = nextVal(p, "o")
      of "min-overlap":
        let v = nextVal(p, "min-overlap")
        try: cfg.minOverlap = parseFloat(v)
        except ValueError:
          stderr.writeLine "error: --min-overlap must be a float, got: " & v
          quit(1)
        cfg.minOverlapSet = true
      of "min-jaccard":
        let v = nextVal(p, "min-jaccard")
        try: cfg.minJaccard = parseFloat(v)
        except ValueError:
          stderr.writeLine "error: --min-jaccard must be a float, got: " & v
          quit(1)
        cfg.minJaccardSet = true
      of "bnd-slop":
        let v = nextVal(p, "bnd-slop")
        try: cfg.bndSlop = parseInt(v)
        except ValueError:
          stderr.writeLine "error: --bnd-slop must be an integer, got: " & v
          quit(1)
        if cfg.bndSlop <= 0:
          stderr.writeLine "error: --bnd-slop must be > 0"
          quit(1)
      of "overwrite":
        cfg.overwrite = true
      of "threads":
        let v = nextVal(p, "threads")
        try: cfg.nThreads = parseInt(v)
        except ValueError:
          stderr.writeLine "error: --threads must be an integer, got: " & v
          quit(1)
        if cfg.nThreads < 1:
          stderr.writeLine "error: --threads must be >= 1"
          quit(1)
      of "tmp-dir":
        cfg.tmpDir = nextVal(p, "tmp-dir")
      of "v", "verbose":
        setVerbose(true)
      of "h", "help":
        annoHelp()
      else:
        stderr.writeLine "error: unknown option: --" & p.key
        annoUsage()
    of cmdArgument:
      positionals.add(p.key)

  if cfg.minOverlapSet == cfg.minJaccardSet:
    if cfg.minOverlapSet:
      stderr.writeLine "error: --min-overlap and --min-jaccard are mutually exclusive"
    else:
      stderr.writeLine "error: exactly one of --min-overlap or --min-jaccard is required"
    annoUsage()
  if positionals.len != 2:
    stderr.writeLine "error: expected 2 input files (input database), got " &
      $positionals.len
    annoUsage()
  cfg.callsetA = positionals[0]
  cfg.callsetB = positionals[1]
  if not fileExists(cfg.callsetA):
    stderr.writeLine "error: input file not found: " & cfg.callsetA
    quit(1)
  if not fileExists(cfg.callsetB):
    stderr.writeLine "error: database file not found: " & cfg.callsetB
    quit(1)
  if cfg.tmpDir == "":
    cfg.tmpDir = getTempDir()

  runAnno(cfg)

proc mainEntry*() =
  var args = commandLineParams()
  while args.len > 0 and (args[0] == "-v" or args[0] == "--verbose"):
    setVerbose(true)
    args.delete(0)
  if args.len == 0:
    usage()
  case args[0]
  of "match":
    runMatch(args[1 .. ^1])
  of "anno":
    runAnnoCli(args[1 .. ^1])
  of "--version":
    echo "matcha v" & VERSION
  of "--help", "-h":
    usage(0)
  else:
    stderr.writeLine "error: unknown subcommand '" & args[0] & "'"
    usage()

when isMainModule:
  mainEntry()
