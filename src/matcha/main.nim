## main.nim — CLI argument parsing and subcommand dispatch.
## Entry point is src/matcha.nim which includes this file.

import std/[os, parseopt, strutils]
import utils, match, log

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
  f.writeLine "Usage: matcha match [options] callsetA callsetB"
  f.writeLine ""
  f.writeLine "Inputs may be VCF.gz (.vcf.gz) or BCF (.bcf); format is detected automatically."
  f.writeLine ""
  f.writeLine "Options:"
  f.writeLine "  --min-overlap FLOAT             minimum reciprocal overlap (0.0-1.0)"
  f.writeLine "  --min-jaccard FLOAT             minimum Jaccard index (0.0-1.0)"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  --output PATH                   output file (default: stdout)"
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      show this help"
  f.writeLine ""
  f.writeLine "At least one of --min-overlap or --min-jaccard is required."
  f.writeLine ""
  f.writeLine "Output is tab-separated with a #-prefixed header line. Columns:"
  f.writeLine "  #CHROM  POS_A  END_A  ID_A  POS_B  END_B  ID_B  SVTYPE  OVERLAP  JACCARD"
  quit(code)

proc runMatch(rawArgs: seq[string]) =
  var cfg = MatchConfig(nThreads: 1)
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
      of "v", "verbose":
        setVerbose(true)
      of "h", "help":
        matchUsage(0)
      else:
        stderr.writeLine "error: unknown option: --" & p.key
        matchUsage()
    of cmdArgument:
      positionals.add(p.key)

  if not cfg.minOverlapSet and not cfg.minJaccardSet:
    stderr.writeLine "error: at least one of --min-overlap or --min-jaccard is required"
    matchUsage()
  if positionals.len != 2:
    stderr.writeLine "error: expected 2 input files, got " & $positionals.len
    matchUsage()

  cfg.callsetA = positionals[0]
  cfg.callsetB = positionals[1]

  if not fileExists(cfg.callsetA):
    stderr.writeLine "error: input file not found: " & cfg.callsetA
    quit(1)
  if not fileExists(cfg.callsetB):
    stderr.writeLine "error: input file not found: " & cfg.callsetB
    quit(1)

  if cfg.tmpDir == "":
    cfg.tmpDir = getTempDir()

  runMatch(cfg)

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
  of "--version":
    echo "matcha v" & VERSION
  of "--help", "-h":
    usage(0)
  else:
    stderr.writeLine "error: unknown subcommand '" & args[0] & "'"
    usage()

when isMainModule:
  mainEntry()
