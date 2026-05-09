## test_utils.nim — shared test helpers: timed template + watchdog thread

import std/[times, atomics, os, locks, strformat, strutils, posix]

when defined(nimPanics):
  {.warning: "tests compiled with --panics:on; AssertionDefect is uncatchable and FAIL line won't print. Set `switch(\"panics\", \"off\")` in config.nims.".}

let gTimeoutSec* = block:
  let e = getEnv("MATCHA_TEST_TIMEOUT", "10")
  try: parseFloat(e) except ValueError: 10.0

let gIgnoreTimeout* = getEnv("MATCHA_IGNORE_TIMEOUT", "0") == "1"

var gDeadline: Atomic[float64]
var gLabel: string
var gLock: Lock
initLock(gLock)

# Current test context — read by onUnhandledException when a test crashes.
var gCurrentId: string
var gCurrentDesc: string
var gCurrentT0: float64

proc watchdog {.thread, gcsafe.} =
  while true:
    sleep(50)
    let dl = gDeadline.load(moRelaxed)
    if dl > 0.0 and epochTime() > dl:
      var msg: string
      withLock(gLock):
        {.gcsafe.}: msg = gLabel & &"\tFAIL\tTIMEOUT (>{gTimeoutSec} s)"
      stderr.writeLine(msg)
      quit(1)

# Only start watchdog when not using fork-based timeout isolation
when true:
  if not gIgnoreTimeout:
    var gWatchdogThread: Thread[void]
    createThread(gWatchdogThread, watchdog)

# Hook invoked on any unhandled exception (including AssertionDefect when
# --panics:off). Runs before Nim's default stack-trace-and-quit behavior.
system.onUnhandledException = proc(errorMsg: string) {.nimcall, gcsafe, raises: [].} =
  try:
    {.gcsafe.}:
      let el = epochTime() - gCurrentT0
      stderr.writeLine(gCurrentId & "\tFAIL\t" &
                        formatFloat(el, ffDecimal, 2) & " s\t" & gCurrentDesc)
      stderr.writeLine(errorMsg.strip)
  except:
    discard

template timed*(id, desc: string, body: untyped) =
  block:
    gCurrentId = id
    gCurrentDesc = desc
    gCurrentT0 = epochTime()

    if gIgnoreTimeout:
      # Fork-based isolation: child runs the body, parent waits with timeout.
      let pid = fork()
      if pid == 0:
        body
        quit(0)
      elif pid > 0:
        var status: cint
        let dl = epochTime() + gTimeoutSec
        var finished = false
        while epochTime() < dl:
          if waitpid(pid, status, WNOHANG) > 0:
            finished = true
            break
          sleep(50)
        if not finished:
          discard kill(pid, cint(SIGKILL))
          discard waitpid(pid, status, 0)
          stderr.writeLine(id & " - " & desc & &"\tFAIL\tTIMEOUT (>{gTimeoutSec} s)")
        elif WIFEXITED(status) and WEXITSTATUS(status) == 0:
          let el = epochTime() - gCurrentT0
          echo id & "\tPASS\t" & formatFloat(el, ffDecimal, 2) & " s\t" & desc
        else:
          let el = epochTime() - gCurrentT0
          stderr.writeLine(id & "\tFAIL\t" &
                            formatFloat(el, ffDecimal, 2) & " s\t" & desc)
      else:
        stderr.writeLine(id & "\tFAIL\tfork() failed")
    else:
      # Direct execution with watchdog thread timeout
      withLock(gLock):
        gLabel = id & " - " & desc
      gDeadline.store(epochTime() + gTimeoutSec, moRelaxed)
      body
      gDeadline.store(0.0, moRelaxed)
      let el = epochTime() - gCurrentT0
      echo id & "\tPASS\t" & formatFloat(el, ffDecimal, 2) & " s\t" & desc
