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
import tables
include "base"


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
  var targets = initTable[string, Target]()

  for t in b.hdr.targets:
    targets[t.name] = t

  let N = max(lc[ x.tid | (x <- b.hdr.targets), int]) + 1
  var M = zeros[int](N, N)
  debug("Initiating matrix of size $# x $#".format(N, N))

  for record in b:
    var
       qi = targets[record.chrom].tid
       mi = targets[record.mate_chrom].tid

    if qi > mi:
      swap(qi, mi)

    M[qi, mi] = M[qi, mi] + 1
    mappings.add((record.qname, record.chrom))

  mappings.sort do (x, y: ReadMapping) -> int:
    result = cmp(x.name, y.name)

  for mapping in mappings[0..<100]:
    echo mapping