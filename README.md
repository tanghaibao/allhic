# ALLHIC Genome scaffolding based on HiC data

[![Travis-CI](https://travis-ci.org/tanghaibao/allhic.svg?branch=master)](https://travis-ci.org/tanghaibao/allhic)

**This software is currently under active development. DO NOT USE.**

## Usage

### Partition

Given a target `k`, number of partitions, the goal of the partitioning
is to separate all the contigs into separate clusters. As with all
clustering algorithm, there is an optimization goal here. The
LACHESIS algorithm is a hierarchical clustering algorithm using
average links.

```bash
allhic partition tests/prunning.sub.bam
```

### Optimize

Given a set of Hi-C contacts between contigs, as specified in the
clmfile, reconstruct the highest scoring ordering and orientations
for these contigs.

```bash
allhic optimize tests/test.clm
```
