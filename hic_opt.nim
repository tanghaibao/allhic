import hic/clm
import strutils
import tables

proc main() {.discardable.} =
  var c = initCLMFile("test", "tests/test.clm")
  c.parse_ids(true)
  c.parse_clm()
  echo c.contacts


when isMainModule:
  main()
