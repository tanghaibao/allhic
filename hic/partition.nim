#
# ALLHIC: Partition function
# Given a target k, number of partitions, the goal of the partitioning is to
# separate all the contigs into separate clusters. As with all clustering
# algorithm, there is an optimization goal here.
#
# The LACHESIS algorithm is a hierarchical clustering algorithm using average
# links.
#


import logging

var logger = newConsoleLogger()
logger.addHandler()

type
  Partitioner* = ref object
    bamfile*: string


proc initPartitioner*(bamfile: string): Partitioner =
  new result
  result.bamfile = bamfile
