#
# ALLHIC: Partition function
# Given a target k, number of partitions, the goal of the partitioning is to
# separate all the contigs into separate clusters. As with all clustering
# algorithm, there is an optimization goal here.
#
# The LACHESIS algorithm is a hierarchical clustering algorithm using average
# links.
#


import algorithm
import hts
import logging

var logger = newConsoleLogger()
logger.addHandler()

type
  Partitioner* = ref object
    bamfile*: string
  ReadMapping* = tuple[name: string, contig: string]


proc initPartitioner*(bamfile: string): Partitioner =
  new result
  result.bamfile = bamfile


proc count_links*(this: Partitioner) =
  var b: Bam
  open(b, this.bamfile, index=false)

  var mappings: seq[ReadMapping] = @[]

  for record in b:
    mappings.add((record.qname, record.chrom))

  mappings.sort do (x, y: ReadMapping) -> int:
    result = cmp(x.name, y.name)

  for mapping in mappings[0..<100]:
    echo mapping
