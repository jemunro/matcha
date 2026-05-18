## log.nim — verbose logging to stderr, gated by a global flag.
##
## All matcha components import this and call logV(). The flag is set once
## at startup by CLI parsing and only read thereafter, so it is safe to read
## from worker threads without synchronisation.

import std/[os, strutils, times]

var gVerbose = false
let gStart = epochTime()

proc setVerbose*(b: bool) =
  gVerbose = b

proc elapsedTag(level: string): string =
  "[" & level & " " & formatFloat(epochTime() - gStart, ffDecimal, 3) & "s] "

proc logV*(msg: string) =
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
