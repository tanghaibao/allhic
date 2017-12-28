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


proc print*[T](x: Matrix[T], filename: string) =
  ## "Pretty prints" the matrix.  All elements of the matrix are
  ## included so large matrices will result in large strings.
  var fw = open(filename, fmWrite)
  debug("Write matrix to `$#`".format(filename))
  for r in countup(0,x.rows()-1):
    for c in countup(0,x.cols()-1):
      fw.write $x[r,c]
      if c != (x.cols()-1):
         fw.write ","
    fw.write "\n"
  fw.close()


proc count_links*(this: Partitioner): Matrix[int] =
  var b: Bam
  open(b, this.bamfile, index=false)

  var targets = b.hdr.targets
  var tig_to_id = initTable[string, int]()
  var id_to_tig = initTable[int, string]()

  for t in targets:
    tig_to_id[t.name] = t.tid
    id_to_tig[t.tid]  = t.name

  let N = max(lc[ x.tid | (x <- b.hdr.targets), int]) + 1
  result = zeros[int](N, N)
  debug("Initiating matrix of size $# x $#".format(N, N))

  for record in b:
    var
       qi = tig_to_id[record.chrom]
       mi = tig_to_id[record.mate_chrom]

    result[qi, mi] = result[qi, mi] + 1
    result[mi, qi] = result[mi, qi] + 1

  result.print("matrix.txt")


proc normalize*(m: Matrix[int]): Matrix[float] =
  discard