#
# ALLHIC: genome scaffolding based on Hi-C data
# Author: Haibao Tang <tanghaibao at gmail dot com>
#

import docopt
import hic/clm
import hic/partition
import strutils
import tables

let doc = """
ALLHIC: genome scaffolding based on Hi-C data

Usage:
  allhic partition
  allhic optimize

Options:
  -h --help       Show this screen.
  --version       Show version.
"""

proc optimize_main() =
  var c = initCLMFile("test", "tests/test.clm")
  c.parse()


proc partition_main() =
  var c = initPartitioner("tests/prunning.sub.bam")
  var M = c.count_links()
  discard M.normalize()


proc main() =
  let args = docopt(doc, version="ALLHIC 0.7.11")

  if args["partition"]:
    partition_main()

  if args["optimize"]:
    optimize_main()


when isMainModule:
  main()
