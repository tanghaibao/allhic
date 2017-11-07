import hic/clm
import strutils
import tables

proc main() {.discardable.} =
  var c = initCLMFile("test", "tests/test.clm")
  c.parse()


when isMainModule:
  main()
