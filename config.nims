# config.nims — automatically loaded by nim for all compilations in this project.
when defined(nimscript):
  import std/[os, strutils]
  # Tests: panics off so AssertionDefect is raised (catchable) instead of
  # triggering rawQuit, letting the `timed` template print a FAIL line.
  if projectName().startsWith("test_"):
    switch("panics", "off")
  # Path to vendored hts-nim (makes `import hts` resolve everywhere).
  switch("path", "vendor/hts-nim/src")
  # Path to src/ (makes `import matcha/utils` etc. resolve from tests/).
  switch("path", "src")
  # Suppress HoleEnumConv warning from hts-nim's Status enum conversions
  # (Status has holes; the conversion is internal to vendored hts-nim and
  # we have no way to patch it from here).
  switch("warning", "HoleEnumConv:off")
# begin Nimble config (version 2)
when withDir(thisDir(), system.fileExists("nimble.paths")):
  include "nimble.paths"
# end Nimble config
