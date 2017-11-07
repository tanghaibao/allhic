import hic/clm
import strutils
import tables

proc main() {.discardable.} =
  var c = initCLMFile("test", "tests/test.clm")
  c.parse_ids(true)
  for key, val in c.tig_to_size.pairs():
    echo("$# $#".format(key, val))
  #c.parse()


when isMainModule:
  main()
