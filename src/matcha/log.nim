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

proc verbose*(): bool =
  gVerbose

proc logV*(msg: string) =
  if not gVerbose: return
  let el = epochTime() - gStart
  stderr.writeLine("[" & formatFloat(el, ffDecimal, 3) & "s] " & msg)

const PreprocWarnPrefix* = "[matcha preproc WARN] "

proc logWarn*(msg: string) =
  ## Emit a warning unconditionally to stderr (independent of -v/--verbose).
  stderr.writeLine(PreprocWarnPrefix & msg)

proc warnCap*(): int =
  ## Per-reason cap for per-record warnings. Override via MATCHA_WARN_CAP.
  let s = getEnv("MATCHA_WARN_CAP", "5")
  try: parseInt(s) except ValueError: 5
