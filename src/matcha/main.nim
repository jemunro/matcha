## main.nim — CLI argument parsing and subcommand dispatch.
## Entry point is src/matcha.nim which includes this file.

import std/[os, parseopt, sequtils, strutils]
import utils, match, anno, collapse, merge, mergecore, preproc, log

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
  logError("--" & flag & " requires a value")
  quit(1)

proc parseFloatOpt(v, flag: string): float64 =
  try: parseFloat(v)
  except ValueError:
    logError("--" & flag & " must be a float, got: " & v); quit(1)

proc parseIntOpt(v, flag: string): int =
  try: parseInt(v)
  except ValueError:
    logError("--" & flag & " must be an integer, got: " & v); quit(1)

proc usage(code: int = 1) =
  let f = if code == 0: stdout else: stderr
  f.writeLine "matcha v" & VERSION
  f.writeLine ""
  f.writeLine "Usage: matcha <subcommand> [options]"
  f.writeLine ""
  f.writeLine "Subcommands:"
  f.writeLine "  match      Find pairwise matches between two SV callsets"
  f.writeLine "  anno       Annotate a callset with INFO fields from a database VCF"
  f.writeLine "  collapse   Cluster equivalent SVs from multiple callers, emit one representative"
  f.writeLine "  merge      Merge per-sample SV callsets into a cohort multi-sample pVCF"
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
  f.writeLine "  --bnd-slop INT                  max breakend offset for BND matches (default: 50)"
  f.writeLine "  --self                          match a single input against itself"
  f.writeLine "                                  (each pair emitted once; no self-self)"
  f.writeLine "  --info FIELDS                   comma-separated INFO fields to add as INFO_A/INFO_B columns"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  --output PATH                   output file (default: stdout)"
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      show this help"
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
  f.writeLine "BND rows use --bnd-slop independently."
  f.writeLine ""
  f.writeLine "Output is tab-separated. A ##matcha_metric=<overlap|jaccard> preamble"
  f.writeLine "precedes the #-prefixed header line. Columns:"
  f.writeLine "  #CHROM  POS_A  END_A  ID_A  POS_B  END_B  ID_B  SVTYPE  SIMILARITY"
  f.writeLine "BND rows emit '.' for END_A / END_B."
  quit(code)

proc runMatch(rawArgs: seq[string]) =
  var cfg = MatchConfig(nThreads: 1, bndSlop: 50, metric: mJaccard, threshold: 0.75)
  var positionals: seq[string]
  # Track which metric flag(s) were supplied so we can enforce the xor rule
  # after parsing. The cfg.metric / cfg.threshold fields don't carry "unset"
  # state on their own.
  var overlapSet, jaccardSet: bool
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      case p.key
      of "min-overlap":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-overlap"), "min-overlap")
        cfg.metric = mOverlap; overlapSet = true
      of "min-jaccard":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-jaccard"), "min-jaccard")
        cfg.metric = mJaccard; jaccardSet = true
      of "bnd-slop":
        cfg.bndSlop = parseIntOpt(nextVal(p, "bnd-slop"), "bnd-slop")
        if cfg.bndSlop <= 0:
          logError("--bnd-slop must be > 0"); quit(1)
      of "threads":
        cfg.nThreads = parseIntOpt(nextVal(p, "threads"), "threads")
        if cfg.nThreads < 1:
          logError("--threads must be >= 1"); quit(1)
      of "tmp-dir":
        cfg.tmpDir = nextVal(p, "tmp-dir")
      of "output":
        cfg.outputPath = nextVal(p, "output")
      of "info":
        cfg.infoFields = nextVal(p, "info").split(',').mapIt(it.strip)
      of "self":
        cfg.selfMode = true
      of "v", "verbose":
        setVerbose(true)
      of "h", "help":
        matchUsage(0)
      else:
        logError("unknown option: --" & p.key)
        matchUsage()
    of cmdArgument:
      positionals.add(p.key)

  if overlapSet and jaccardSet:
    logError("--min-overlap and --min-jaccard are mutually exclusive")
    matchUsage()
  let expected = if cfg.selfMode: 1 else: 2
  if positionals.len != expected:
    let what = if cfg.selfMode: "1 input file (--self mode)" else: "2 input files"
    logError("expected " & what & ", got " & $positionals.len)
    matchUsage()

  cfg.callsetA = positionals[0]
  if not cfg.selfMode:
    cfg.callsetB = positionals[1]

  if not fileExists(cfg.callsetA):
    logError("input file not found: " & cfg.callsetA)
    quit(1)
  if not cfg.selfMode and not fileExists(cfg.callsetB):
    logError("input file not found: " & cfg.callsetB)
    quit(1)

  if cfg.tmpDir == "":
    cfg.tmpDir = getTempDir()
  cfg.tmpDir = makeRunTmpDir(cfg.tmpDir)

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
  f.writeLine "  --bnd-slop INT                  max breakend offset for BND matches (default: 50)"
  f.writeLine "  --overwrite                     replace OUTFIELDs already in input header"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      show full help (functions, MATCHA_* variables)"
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
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
  f.writeLine "  --bnd-slop INT                  max breakend offset for BND matches (default: 50)"
  f.writeLine "  --overwrite                     replace OUTFIELDs that already exist in input header"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      show this help"
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
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
  var cfg = AnnoConfig(nThreads: 1, bndSlop: 50, metric: mJaccard, threshold: 0.75)
  var positionals: seq[string]
  var overlapSet, jaccardSet: bool
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
          logError(e.msg)
          quit(1)
      of "o", "output":
        cfg.outputPath = nextVal(p, "o")
      of "min-overlap":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-overlap"), "min-overlap")
        cfg.metric = mOverlap; overlapSet = true
      of "min-jaccard":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-jaccard"), "min-jaccard")
        cfg.metric = mJaccard; jaccardSet = true
      of "bnd-slop":
        cfg.bndSlop = parseIntOpt(nextVal(p, "bnd-slop"), "bnd-slop")
        if cfg.bndSlop <= 0:
          logError("--bnd-slop must be > 0"); quit(1)
      of "overwrite":
        cfg.overwrite = true
      of "threads":
        cfg.nThreads = parseIntOpt(nextVal(p, "threads"), "threads")
        if cfg.nThreads < 1:
          logError("--threads must be >= 1"); quit(1)
      of "tmp-dir":
        cfg.tmpDir = nextVal(p, "tmp-dir")
      of "v", "verbose":
        setVerbose(true)
      of "h", "help":
        annoHelp()
      else:
        logError("unknown option: --" & p.key)
        annoUsage()
    of cmdArgument:
      positionals.add(p.key)

  if overlapSet and jaccardSet:
    logError("--min-overlap and --min-jaccard are mutually exclusive")
    annoUsage()
  if positionals.len != 2:
    logError("expected 2 input files (input database), got " & $positionals.len)
    annoUsage()
  cfg.callsetA = positionals[0]
  cfg.callsetB = positionals[1]
  if not fileExists(cfg.callsetA):
    logError("input file not found: " & cfg.callsetA)
    quit(1)
  if not fileExists(cfg.callsetB):
    logError("database file not found: " & cfg.callsetB)
    quit(1)
  if cfg.tmpDir == "":
    cfg.tmpDir = getTempDir()
  cfg.tmpDir = makeRunTmpDir(cfg.tmpDir)

  runAnno(cfg)

proc collapseUsage(code: int = 1) =
  let f = if code == 0: stdout else: stderr
  f.writeLine "Usage:"
  f.writeLine "  matcha collapse [options] [Name:]callset1.bcf [Name:]callset2.bcf ..."
  f.writeLine ""
  f.writeLine "Cluster equivalent SVs from N single-sample callsets and emit one"
  f.writeLine "representative record per cluster with provenance INFO fields."
  f.writeLine ""
  f.writeLine "Options:"
  f.writeLine "  --min-overlap FLOAT           minimum reciprocal overlap (0.0-1.0)"
  f.writeLine "  --min-jaccard FLOAT           minimum Jaccard index (0.0-1.0)"
  f.writeLine "  --bnd-slop INT                max breakend offset for BND (default: 50)"
  f.writeLine "  --linkage average|single|complete  agglomerative linkage (default: average)"
  f.writeLine "  --priority CRITERIA           comma-separated: PASS,QUAL,CENTRE,ORDER"
  f.writeLine "                                default: PASS,CENTRE,ORDER"
  f.writeLine "                                ORDER is always appended as final tiebreaker"
  f.writeLine "  --format FIELDS               comma-separated FORMAT fields to carry"
  f.writeLine "                                default: GT"
  f.writeLine "  --info FIELDS                 comma-separated INFO fields to keep"
  f.writeLine "                                default: SVTYPE,SVLEN,END,CHR2,POS2 only"
  f.writeLine "  -o, --output PATH             output file (.vcf | .vcf.gz | .bcf)"
  f.writeLine "                                default: uncompressed VCF to stdout"
  f.writeLine "  --threads INT                 worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                temp directory (default: system temp)"
  f.writeLine "  -v, --verbose                 verbose logging to stderr"
  f.writeLine "  -h, --help                    show this help"
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
  f.writeLine ""
  f.writeLine "Input names: positional args may be prefixed with 'Name:' (e.g. Delly:delly.bcf)."
  f.writeLine "If no name prefix is given, the basename without extension is used."
  f.writeLine ""
  f.writeLine "Output INFO fields added: CALLERS, N_CALLERS, N_MERGED."
  quit(code)

proc parsePriority(s: string): seq[PriorityCriterion] =
  for tok in s.split(','):
    case tok.strip.toUpperAscii
    of "PASS":   result.add(pcPass)
    of "QUAL":   result.add(pcQual)
    of "CENTRE", "CENTER": result.add(pcCentre)
    of "ORDER":  result.add(pcOrder)
    else:
      logError("unknown priority criterion '" & tok & "' — valid values: PASS, QUAL, CENTRE, ORDER")
      quit(1)
  # ORDER is always the final tiebreaker; append if not already last.
  if result.len == 0 or result[^1] != pcOrder:
    result.add(pcOrder)

proc runCollapseCli(rawArgs: seq[string]) =
  var cfg = CollapseConfig(
    nThreads:     1,
    bndSlop:      50,
    metric:       mJaccard,
    threshold:    0.75,
    linkage:      lmAverage,
    priority:     @[pcPass, pcCentre, pcOrder],
    formatFields: @["GT"],
  )
  var positionals: seq[string]
  var overlapSet, jaccardSet: bool
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      case p.key
      of "min-overlap":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-overlap"), "min-overlap")
        cfg.metric = mOverlap; overlapSet = true
      of "min-jaccard":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-jaccard"), "min-jaccard")
        cfg.metric = mJaccard; jaccardSet = true
      of "bnd-slop":
        cfg.bndSlop = parseIntOpt(nextVal(p, "bnd-slop"), "bnd-slop")
        if cfg.bndSlop <= 0:
          logError("--bnd-slop must be > 0"); quit(1)
      of "linkage":
        let v = nextVal(p, "linkage").toLowerAscii
        case v
        of "average":  cfg.linkage = lmAverage
        of "single":   cfg.linkage = lmSingle
        of "complete": cfg.linkage = lmComplete
        else:
          logError("--linkage must be average, single, or complete")
          quit(1)
      of "priority":
        cfg.priority = parsePriority(nextVal(p, "priority"))
      of "format":
        cfg.formatFields = nextVal(p, "format").split(',').mapIt(it.strip)
      of "info":
        cfg.infoFields = nextVal(p, "info").split(',').mapIt(it.strip)
      of "o", "output":
        cfg.outputPath = nextVal(p, "o")
      of "threads":
        cfg.nThreads = parseIntOpt(nextVal(p, "threads"), "threads")
        if cfg.nThreads < 1:
          logError("--threads must be >= 1"); quit(1)
      of "tmp-dir":
        cfg.tmpDir = nextVal(p, "tmp-dir")
      of "v", "verbose": setVerbose(true)
      of "h", "help":    collapseUsage(0)
      else:
        logError("unknown option: --" & p.key)
        collapseUsage()
    of cmdArgument:
      positionals.add(p.key)

  if overlapSet and jaccardSet:
    logError("--min-overlap and --min-jaccard are mutually exclusive")
    collapseUsage()

  if positionals.len < 1:
    logError("at least one input file is required")
    collapseUsage()

  # Parse [Name:]path positional arguments.
  for arg in positionals:
    let colonIdx = arg.find(':')
    var name, path: string
    if colonIdx >= 0:
      name = arg[0 ..< colonIdx]
      path = arg[colonIdx + 1 .. ^1]
    else:
      path = arg
      name = splitFile(path).name
    if name.len == 0: name = splitFile(path).name
    cfg.callers.add(CallerInput(name: name, path: path))

  for caller in cfg.callers:
    if not fileExists(caller.path):
      logError("input file not found: " & caller.path)
      quit(1)

  if cfg.tmpDir == "":
    cfg.tmpDir = getTempDir()
  cfg.tmpDir = makeRunTmpDir(cfg.tmpDir)

  # Build command line string for provenance header.
  let cmdLine = "matcha collapse " & rawArgs.join(" ")

  runCollapse(cfg, cmdLine)

proc mergeUsage(code: int = 1) =
  let f = if code == 0: stdout else: stderr
  f.writeLine "Usage:"
  f.writeLine "  matcha merge [options] [Name:]callset1.bcf [Name:]callset2.bcf ..."
  f.writeLine ""
  f.writeLine "Merge per-sample SV callsets (≥2 inputs, exactly 1 sample each, distinct"
  f.writeLine "sample IDs) into a multi-sample cohort pVCF. One row per cluster; FORMAT"
  f.writeLine "columns per sample; cohort INFO AC/AN/AF computed from assembled GTs."
  f.writeLine ""
  f.writeLine "Options:"
  f.writeLine "  --min-overlap FLOAT           minimum reciprocal overlap (0.0-1.0)"
  f.writeLine "  --min-jaccard FLOAT           minimum Jaccard index (0.0-1.0)"
  f.writeLine "  --bnd-slop INT                max breakend offset for BND (default: 50)"
  f.writeLine "  --linkage average|single|complete  agglomerative linkage (default: average)"
  f.writeLine "  --priority CRITERIA           comma-separated: PASS,QUAL,CENTRE,ORDER"
  f.writeLine "                                default: PASS,CENTRE,ORDER (drives representative)"
  f.writeLine "  --format FIELDS               comma-separated FORMAT fields to carry"
  f.writeLine "                                default: GT (auto-added if absent)"
  f.writeLine "  --info FIELDS                 comma-separated INFO fields to keep from rep"
  f.writeLine "                                default: only auto-extracted + cohort + CALLERS"
  f.writeLine "  -o, --output PATH             output file (.vcf | .vcf.gz | .bcf)"
  f.writeLine "                                default: uncompressed VCF to stdout"
  f.writeLine "  --threads INT                 worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                temp directory (default: system temp)"
  f.writeLine "  -v, --verbose                 verbose logging to stderr"
  f.writeLine "  -h, --help                    show this help"
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
  f.writeLine ""
  f.writeLine "Input names: positional args may be prefixed with 'Name:' (e.g. S1:s1.bcf)."
  f.writeLine "If no name prefix is given, the basename without extension is used."
  f.writeLine ""
  f.writeLine "Output INFO fields added: AC, AN, AF (always); CALLERS, N_CALLERS (when inputs had them)."
  quit(code)

proc runMergeCli(rawArgs: seq[string]) =
  var cfg = MergeConfig(
    nThreads:     1,
    bndSlop:      50,
    metric:       mJaccard,
    threshold:    0.75,
    linkage:      lmAverage,
    priority:     @[pcPass, pcCentre, pcOrder],
    formatFields: @["GT"],
  )
  var positionals: seq[string]
  var overlapSet, jaccardSet: bool
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      case p.key
      of "min-overlap":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-overlap"), "min-overlap")
        cfg.metric = mOverlap; overlapSet = true
      of "min-jaccard":
        cfg.threshold = parseFloatOpt(nextVal(p, "min-jaccard"), "min-jaccard")
        cfg.metric = mJaccard; jaccardSet = true
      of "bnd-slop":
        cfg.bndSlop = parseIntOpt(nextVal(p, "bnd-slop"), "bnd-slop")
        if cfg.bndSlop <= 0:
          logError("--bnd-slop must be > 0"); quit(1)
      of "linkage":
        let v = nextVal(p, "linkage").toLowerAscii
        case v
        of "average":  cfg.linkage = lmAverage
        of "single":   cfg.linkage = lmSingle
        of "complete": cfg.linkage = lmComplete
        else:
          logError("--linkage must be average, single, or complete")
          quit(1)
      of "priority":
        cfg.priority = parsePriority(nextVal(p, "priority"))
      of "format":
        cfg.formatFields = nextVal(p, "format").split(',').mapIt(it.strip)
      of "info":
        cfg.infoFields = nextVal(p, "info").split(',').mapIt(it.strip)
      of "o", "output":
        cfg.outputPath = nextVal(p, "o")
      of "threads":
        cfg.nThreads = parseIntOpt(nextVal(p, "threads"), "threads")
        if cfg.nThreads < 1:
          logError("--threads must be >= 1"); quit(1)
      of "tmp-dir":
        cfg.tmpDir = nextVal(p, "tmp-dir")
      of "v", "verbose": setVerbose(true)
      of "h", "help":    mergeUsage(0)
      else:
        logError("unknown option: --" & p.key)
        mergeUsage()
    of cmdArgument:
      positionals.add(p.key)

  if overlapSet and jaccardSet:
    logError("--min-overlap and --min-jaccard are mutually exclusive")
    mergeUsage()
  if positionals.len < 2:
    logError("merge requires at least 2 input files (got " & $positionals.len & ")")
    mergeUsage()

  for arg in positionals:
    let colonIdx = arg.find(':')
    var name, path: string
    if colonIdx >= 0:
      name = arg[0 ..< colonIdx]
      path = arg[colonIdx + 1 .. ^1]
    else:
      path = arg
      name = splitFile(path).name
    if name.len == 0: name = splitFile(path).name
    cfg.callers.add(CallerInput(name: name, path: path))

  for caller in cfg.callers:
    if not fileExists(caller.path):
      logError("input file not found: " & caller.path)
      quit(1)

  if cfg.tmpDir == "":
    cfg.tmpDir = getTempDir()
  cfg.tmpDir = makeRunTmpDir(cfg.tmpDir)

  let cmdLine = "matcha merge " & rawArgs.join(" ")
  runMerge(cfg, cmdLine)

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
  of "collapse":
    runCollapseCli(args[1 .. ^1])
  of "merge":
    runMergeCli(args[1 .. ^1])
  of "--version":
    echo "matcha v" & VERSION
  of "--help", "-h":
    usage(0)
  else:
    logError("unknown subcommand '" & args[0] & "'")
    usage()

when isMainModule:
  mainEntry()
