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


proc initPartitioner*(bamfile: string): Partitioner =
  new result
  result.bamfile = bamfile


proc count_links*(this: Partitioner) =
  var b: Bam
  open(b, this.bamfile, index=false)

  var targets = b.hdr.targets
  var tig_to_id = initTable[string, int]()
  var id_to_tig = initTable[int, string]()

  for t in targets:
    tig_to_id[t.name] = t.tid
    id_to_tig[t.tid]  = t.name

  let N = max(lc[ x.tid | (x <- b.hdr.targets), int]) + 1
  var M = zeros[int](N, N)
  debug("Initiating matrix of size $# x $#".format(N, N))

  for record in b:
    var
       qi = tig_to_id[record.chrom]
       mi = tig_to_id[record.mate_chrom]

    M[qi, mi] = M[qi, mi] + 1
    M[mi, qi] = M[mi, qi] + 1

  for i in 0..<N:
    var itig = id_to_tig[i]
    for j in i..<N:
      var jtig = id_to_tig[j]
      if M[i, j] > 1:
        echo format("M[$#, $#] = $#", itig, jtig, M[i, j])
