## main.nim — CLI argument parsing and subcommand dispatch.
## Entry point is src/matcha.nim which includes this file.

import std/[os, parseopt, sequtils, sets, strutils]
import utils, match, anno, collapse, merge, mergecore, preproc, log

const VERSION = MatchaVersion

# Embedded at compile time so the release binary is self-contained for
# distribution. Paths are relative to this source file.
const matchaLicense = staticRead("../../LICENSE")
const htsNimLicense = staticRead("../../vendor/hts-nim/LICENSE")

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

proc parseChunkSizeArg*(raw: string): int64 =
  if raw.len == 0:
    logError("--chunk-size requires a value"); quit(1)
  let suffix = raw[^1]
  let digits = raw[0 ..< raw.len - 1]
  if suffix != 'K' and suffix != 'M':
    logError("--chunk-size suffix must be K or M (e.g. 100K, 50M), got: " & raw); quit(1)
  var n: int64
  try: n = parseInt(digits)
  except ValueError:
    logError("--chunk-size prefix must be a positive integer, got: " & raw); quit(1)
  if n <= 0:
    logError("--chunk-size must be > 0, got: " & raw); quit(1)
  result = if suffix == 'K': n * 1_000 else: n * 1_000_000

# Shared CLI state for the flags that are identical across all subcommands.
type SharedOpts = object
  metric:    Metric
  threshold: float64
  bndSlop, insSlop, nThreads: int
  insMinSim: float64
  tmpDir:    string
  keptChrs:  seq[string]   ## --chrs: active set; empty = all input contigs.
  chrSet:    seq[string]   ## --chr-set: universe for BND mates + header contigs;
                           ## empty = all input contigs (no BND-mate filter).
  chunkSize: int64         ## --chunk-size: A-side POS range per job.
  overlapSet, jaccardSet, useShm: bool

proc initSharedOpts(): SharedOpts =
  SharedOpts(metric: mJaccard, threshold: 0.75,
             bndSlop: 50, insSlop: 50, insMinSim: 0.75, nThreads: 1,
             chunkSize: 50_000_000)

proc parseSharedOpt(s: var SharedOpts, p: var OptParser, key: string): bool =
  result = true
  case key
  of "min-overlap":
    s.threshold = parseFloatOpt(nextVal(p, "min-overlap"), "min-overlap")
    if s.threshold <= 0 or s.threshold > 1:
      logError("--min-overlap must be in (0, 1]"); quit(1)
    s.metric = mOverlap; s.overlapSet = true
  of "min-jaccard":
    s.threshold = parseFloatOpt(nextVal(p, "min-jaccard"), "min-jaccard")
    if s.threshold <= 0 or s.threshold > 1:
      logError("--min-jaccard must be in (0, 1]"); quit(1)
    s.metric = mJaccard; s.jaccardSet = true
  of "bnd-slop":
    s.bndSlop = parseIntOpt(nextVal(p, "bnd-slop"), "bnd-slop")
    if s.bndSlop <= 0:
      logError("--bnd-slop must be > 0"); quit(1)
  of "ins-slop":
    s.insSlop = parseIntOpt(nextVal(p, "ins-slop"), "ins-slop")
    if s.insSlop <= 0:
      logError("--ins-slop must be > 0"); quit(1)
  of "min-ins-sim":
    s.insMinSim = parseFloatOpt(nextVal(p, "min-ins-sim"), "min-ins-sim")
    if s.insMinSim <= 0 or s.insMinSim > 1:
      logError("--min-ins-sim must be in (0, 1]"); quit(1)
  of "threads":
    s.nThreads = parseIntOpt(nextVal(p, "threads"), "threads")
    if s.nThreads < 1:
      logError("--threads must be >= 1"); quit(1)
  of "tmp-dir":
    s.tmpDir = nextVal(p, "tmp-dir")
  of "use-shm":
    s.useShm = true
  of "chrs":
    s.keptChrs = parseChrsArg(nextVal(p, "chrs"))
  of "chr-set":
    s.chrSet = parseChrsArg(nextVal(p, "chr-set"))
  of "chunk-size":
    s.chunkSize = parseChunkSizeArg(nextVal(p, "chunk-size"))
  of "v", "verbose":
    setVerbose(true)
  else:
    result = false

proc resolveTmpDir(explicit: string, useShm: bool): string =
  if explicit != "":
    if useShm: logWarn("--use-shm ignored because --tmp-dir was also given")
    makeRunTmpDir(explicit)
  elif useShm:
    makeShmRunTmpDir()
  else:
    makeRunTmpDir(getTempDir())

template applySharedOpts(s: SharedOpts, cfg: typed) =
  cfg.metric    = s.metric
  cfg.threshold = s.threshold
  cfg.bndSlop   = s.bndSlop
  cfg.insSlop   = s.insSlop
  cfg.insMinSim = s.insMinSim
  cfg.nThreads  = s.nThreads
  cfg.keptChrs  = s.keptChrs
  cfg.chrSet    = s.chrSet
  cfg.chunkSize = s.chunkSize
  # --chr-set given without --chrs: default the active set to chr-set so
  # output records and header contigs stay consistent.
  if cfg.chrSet.len > 0 and cfg.keptChrs.len == 0:
    cfg.keptChrs = cfg.chrSet
  # --chrs must be ⊆ --chr-set when both are provided.
  if cfg.chrSet.len > 0 and s.keptChrs.len > 0:
    let chrSetH = toHashSet(cfg.chrSet)
    for c in s.keptChrs:
      if c notin chrSetH:
        logError("--chrs entry '" & c & "' is not in --chr-set"); quit(1)
  cfg.tmpDir    = resolveTmpDir(s.tmpDir, s.useShm)

# Cluster CLI state for the flags shared between collapse and merge.
type ClusterOpts = object
  linkage:      LinkageMethod
  priority:     seq[PriorityCriterion]
  formatFields: seq[string]
  infoFields:   seq[string]
  outputPath:   string

proc initClusterOpts(): ClusterOpts =
  ClusterOpts(linkage: lmAverage,
              priority: @[pcPass, pcCentre, pcOrder],
              formatFields: @["GT"])

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

proc parseClusterOpt(c: var ClusterOpts, p: var OptParser, key: string): bool =
  result = true
  case key
  of "linkage":
    let v = nextVal(p, "linkage").toLowerAscii
    case v
    of "average":  c.linkage = lmAverage
    of "single":   c.linkage = lmSingle
    of "complete": c.linkage = lmComplete
    else:
      logError("--linkage must be average, single, or complete")
      quit(1)
  of "priority":
    c.priority = parsePriority(nextVal(p, "priority"))
  of "format":
    c.formatFields = nextVal(p, "format").split(',').mapIt(it.strip)
  of "info":
    c.infoFields = nextVal(p, "info").split(',').mapIt(it.strip)
  of "o", "output":
    c.outputPath = nextVal(p, "o")
  else:
    result = false

template applyClusterOpts(c: ClusterOpts, cfg: typed) =
  cfg.linkage      = c.linkage
  cfg.priority     = c.priority
  cfg.formatFields = c.formatFields
  cfg.infoFields   = c.infoFields
  cfg.outputPath   = c.outputPath

# Column width 34: option text padded to 34 chars before the description.
proc writeMetricUsageLines(f: File) =
  f.writeLine "  --min-overlap FLOAT             minimum reciprocal overlap (0.0-1.0)"
  f.writeLine "  --min-jaccard FLOAT             minimum Jaccard index (0.0-1.0)"
  f.writeLine "  --bnd-slop INT                  max breakend offset for BND (default: 50)"
  f.writeLine "  --min-ins-sim FLOAT             minimum INS combined sim = sqrt(pos*len) (default: 0.75)"
  f.writeLine "  --ins-slop INT                  max position offset for INS (default: 50)"

proc writePreprocUsageLines(f: File) =
  f.writeLine "  --chrs CHR[,CHR...]             restrict output to listed chromosomes (active set)"
  f.writeLine "  --chr-set CHR[,CHR...]          universe: BND-mate filter + header contigs"
  f.writeLine "                                  (default: all input contigs; no BND drop)"
  f.writeLine "  --chunk-size INT[K|M]           A-side POS window per parallel job (default: 50M)"
  f.writeLine "                                  suffix K=×1,000 or M=×1,000,000 required"
  f.writeLine "  --threads INT                   number of worker threads (default: 1)"
  f.writeLine "  --tmp-dir PATH                  temp directory (default: system temp)"
  f.writeLine "  --use-shm                       write temp BCFs to /dev/shm (RAM); faster if $TMPDIR"
  f.writeLine "                                  is on a network filesystem; risks OOM for large inputs"

proc writeVerboseHelpLines(f: File, helpDesc = "show this help") =
  f.writeLine "  -v, --verbose                   verbose logging to stderr"
  f.writeLine "  -h, --help                      " & helpDesc

proc writeClusterUsageLines(f: File,
                            priorityDefault: string,
                            includeOrderNote: bool,
                            formatDefault: string,
                            infoFirstLine: string,
                            infoDefault: string) =
  f.writeLine "  --linkage MODE                  agglomerative clustering linkage (default: average)"
  f.writeLine "                                  one of: average, single, complete"
  f.writeLine "  --priority CRITERIA             comma-separated: PASS,QUAL,CENTRE,ORDER"
  f.writeLine "                                  default: " & priorityDefault
  if includeOrderNote:
    f.writeLine "                                  ORDER is always appended as final tiebreaker"
  f.writeLine "  --format FIELDS                 comma-separated FORMAT fields to carry"
  f.writeLine "                                  default: " & formatDefault
  f.writeLine "  --info FIELDS                   " & infoFirstLine
  f.writeLine "                                  default: " & infoDefault
  f.writeLine "  -o, --output PATH               output file (.vcf | .vcf.gz | .bcf)"
  f.writeLine "                                  default: uncompressed VCF to stdout"

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
  f.writeLine "  --license       show license information"
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
  writeMetricUsageLines(f)
  f.writeLine "  --self                          match a single input against itself"
  f.writeLine "                                  (each pair emitted once; no self-self)"
  f.writeLine "  --info FIELDS                   comma-separated INFO fields to add as INFO_A/INFO_B columns"
  writePreprocUsageLines(f)
  f.writeLine "  --output PATH                   output file (default: stdout)"
  writeVerboseHelpLines(f)
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
  f.writeLine "BND rows use --bnd-slop independently."
  f.writeLine "INS rows use --min-ins-sim with --ins-slop (combined position+size similarity)."
  f.writeLine ""
  f.writeLine "Output is tab-separated. A ##matcha_metric=<overlap|jaccard> preamble"
  f.writeLine "precedes the #-prefixed header line. Columns:"
  f.writeLine "  #CHROM  POS_A  END_A  ID_A  POS_B  END_B  ID_B  SVTYPE  SIMILARITY"
  f.writeLine "BND rows emit '.' for END_A / END_B."
  quit(code)

proc runMatch(rawArgs: seq[string]) =
  if rawArgs.len == 0: matchUsage(0)
  var cfg = MatchConfig()
  var s = initSharedOpts()
  var positionals: seq[string]
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      if parseSharedOpt(s, p, p.key): continue
      case p.key
      of "output":
        cfg.outputPath = nextVal(p, "output")
      of "info":
        cfg.infoFields = nextVal(p, "info").split(',').mapIt(it.strip)
      of "self":
        cfg.selfMode = true
      of "h", "help":
        matchUsage(0)
      else:
        logError("unknown option: --" & p.key)
        matchUsage()
    of cmdArgument:
      positionals.add(p.key)

  if s.overlapSet and s.jaccardSet:
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

  applySharedOpts(s, cfg)
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
  writeMetricUsageLines(f)
  f.writeLine "  --overwrite                     replace OUTFIELDs already in input header"
  writePreprocUsageLines(f)
  writeVerboseHelpLines(f, "show full help (functions, MATCHA_* variables)")
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
  writeMetricUsageLines(f)
  f.writeLine "  --overwrite                     replace OUTFIELDs that already exist in input header"
  writePreprocUsageLines(f)
  writeVerboseHelpLines(f)
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
  if rawArgs.len == 0: annoUsage(0)
  var cfg = AnnoConfig()
  var s = initSharedOpts()
  var positionals: seq[string]
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      if parseSharedOpt(s, p, p.key): continue
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
      of "overwrite":
        cfg.overwrite = true
      of "h", "help":
        annoHelp()
      else:
        logError("unknown option: --" & p.key)
        annoUsage()
    of cmdArgument:
      positionals.add(p.key)

  if s.overlapSet and s.jaccardSet:
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
  applySharedOpts(s, cfg)
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
  writeMetricUsageLines(f)
  writeClusterUsageLines(f,
    priorityDefault  = "PASS,CENTRE,ORDER",
    includeOrderNote = true,
    formatDefault    = "GT",
    infoFirstLine    = "comma-separated INFO fields to keep",
    infoDefault      = "SVTYPE,SVLEN,END,CHR2,POS2 only")
  writePreprocUsageLines(f)
  writeVerboseHelpLines(f)
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
  f.writeLine ""
  f.writeLine "Input names: positional args may be prefixed with 'Name:' (e.g. Delly:delly.bcf)."
  f.writeLine "If no name prefix is given, the basename without extension is used."
  f.writeLine ""
  f.writeLine "Output INFO fields added: CALLERS, N_CALLERS, N_MERGED."
  quit(code)

proc runCollapseCli(rawArgs: seq[string]) =
  if rawArgs.len == 0: collapseUsage(0)
  var cfg = CollapseConfig()
  var s = initSharedOpts()
  var c = initClusterOpts()
  var positionals: seq[string]
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      if parseSharedOpt(s, p, p.key): continue
      if parseClusterOpt(c, p, p.key): continue
      case p.key
      of "h", "help":
        collapseUsage(0)
      else:
        logError("unknown option: --" & p.key)
        collapseUsage()
    of cmdArgument:
      positionals.add(p.key)

  if s.overlapSet and s.jaccardSet:
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

  applySharedOpts(s, cfg)
  applyClusterOpts(c, cfg)

  # Build command line string for provenance header.
  let cmdLine = "matcha collapse " & rawArgs.join(" ")

  runCollapse(cfg, cmdLine)

proc mergeUsage(code: int = 1) =
  let f = if code == 0: stdout else: stderr
  f.writeLine "Usage:"
  f.writeLine "  matcha merge [options] callset1.bcf callset2.bcf ..."
  f.writeLine ""
  f.writeLine "Merge per-sample SV callsets (≥2 inputs, exactly 1 sample each, distinct"
  f.writeLine "sample IDs) into a multi-sample cohort pVCF. One row per cluster; FORMAT"
  f.writeLine "columns per sample; cohort INFO AC/AN/AF computed from assembled GTs."
  f.writeLine ""
  f.writeLine "Options:"
  writeMetricUsageLines(f)
  writeClusterUsageLines(f,
    priorityDefault  = "PASS,CENTRE,ORDER (drives representative)",
    includeOrderNote = false,
    formatDefault    = "GT (auto-added if absent)",
    infoFirstLine    = "comma-separated INFO fields to keep from rep",
    infoDefault      = "only auto-extracted + cohort + CALLERS")
  f.writeLine "  --missing-to-ref                treat absent samples as 0/0 (count toward AN; like bcftools merge)"
  writePreprocUsageLines(f)
  writeVerboseHelpLines(f)
  f.writeLine ""
  f.writeLine "Default metric: --min-jaccard 0.75. Specify --min-overlap or --min-jaccard to override."
  f.writeLine ""
  f.writeLine "Output INFO fields added: AC, AN, AF (always); CALLERS, N_CALLERS (when inputs had them)."
  f.writeLine "AC/AN reflect --missing-to-ref: absent samples are counted as 0/0 when the flag is set."
  quit(code)

proc runMergeCli(rawArgs: seq[string]) =
  if rawArgs.len == 0: mergeUsage(0)
  var cfg = MergeConfig()
  var s = initSharedOpts()
  var c = initClusterOpts()
  var positionals: seq[string]
  var p = initOptParser(rawArgs, shortNoVal = ShortNoVal)
  while true:
    p.next()
    case p.kind
    of cmdEnd: break
    of cmdShortOption, cmdLongOption:
      if parseSharedOpt(s, p, p.key): continue
      if parseClusterOpt(c, p, p.key): continue
      case p.key
      of "missing-to-ref":
        cfg.missingToRef = true
      of "h", "help":
        mergeUsage(0)
      else:
        logError("unknown option: --" & p.key)
        mergeUsage()
    of cmdArgument:
      positionals.add(p.key)

  if s.overlapSet and s.jaccardSet:
    logError("--min-overlap and --min-jaccard are mutually exclusive")
    mergeUsage()
  if positionals.len < 2:
    logError("merge requires at least 2 input files (got " & $positionals.len & ")")
    mergeUsage()

  for arg in positionals:
    if ':' in arg:
      logError("matcha merge does not support 'Name:' prefix on inputs " &
               "(sample IDs come from the input BCF headers); got '" & arg & "'")
      quit(1)
    cfg.callers.add(CallerInput(name: splitFile(arg).name, path: arg))

  for caller in cfg.callers:
    if not fileExists(caller.path):
      logError("input file not found: " & caller.path)
      quit(1)

  applySharedOpts(s, cfg)
  applyClusterOpts(c, cfg)
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
  of "--license":
    echo "=== matcha ==="
    echo ""
    echo matchaLicense
    echo "=== hts-nim ==="
    echo ""
    echo htsNimLicense
  of "--help", "-h":
    usage(0)
  else:
    logError("unknown subcommand '" & args[0] & "'")
    usage()

when isMainModule:
  mainEntry()
