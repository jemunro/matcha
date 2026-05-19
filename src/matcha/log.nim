## log.nim — structured logging to stderr.
##
## logInfo    — always emitted; meaningful progress/result lines a user always wants.
## logVerbose — gated by -v/--verbose; per-file and per-job detail.
## logWarn    — always emitted (suppressible in tests via setQuiet).
## logError   — always emitted; fatal condition messages.
##
## The verbose flag is set once at startup and only read thereafter, so it is
## safe to read from worker threads without synchronisation.

import std/[os, strutils, times]

var gVerbose = false
let gStart = epochTime()

proc setVerbose*(b: bool) =
  gVerbose = b

proc elapsedTag(level: string): string =
  "[" & level & " " & formatFloat(epochTime() - gStart, ffDecimal, 3) & "s] "

proc logInfo*(msg: string) =
  stderr.writeLine(elapsedTag("INFO") & msg)

proc logVerbose*(msg: string) =
  if not gVerbose: return
  stderr.writeLine(elapsedTag("INFO") & msg)

var gQuiet = false

proc setQuiet*(b: bool) =
  ## When set, logWarn becomes a no-op. Used by tests to keep PASS/FAIL
  ## output uncluttered. Production code never calls this.
  gQuiet = b

proc logWarn*(msg: string) =
  ## Emit a warning to stderr (independent of -v/--verbose). Suppressed
  ## when setQuiet(true) has been called (test runs).
  if gQuiet: return
  stderr.writeLine(elapsedTag("WARN") & msg)

proc logError*(msg: string) =
  stderr.writeLine(elapsedTag("ERROR") & msg)

proc warnCap*(): int =
  ## Per-reason cap for per-record warnings. Override via MATCHA_WARN_CAP.
  let s = getEnv("MATCHA_WARN_CAP", "5")
  try: parseInt(s) except ValueError: 5
