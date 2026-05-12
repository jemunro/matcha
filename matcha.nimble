# Package
version     = "0.1.0"
author      = "Jacob Munro"
description = "Structural variant matching and annotation"
license     = "MIT"
srcDir      = "src"
bin         = @["matcha"]

# Dependencies
requires "nim >= 2.0.0"

# ---------------------------------------------------------------------------
# Tasks
# ---------------------------------------------------------------------------

task release, "Build release binary":
  --define:release
  --define:strip
  --panics:on
  setCommand "build"

task test, "Run all tests":
  exec "nimble build -y"
  exec "rm -rf nimcache/tests"
  exec "nim c --hints:off -r tests/test_intervals.nim"
  exec "nim c --hints:off -r tests/test_bins.nim"
  exec "nim c --hints:off -r tests/test_preproc.nim"
  exec "nim c --hints:off -r tests/test_match.nim"
  exec "nim c --hints:off -r tests/test_anno.nim"
  exec "nim c --hints:off -r tests/test_collapse.nim"
