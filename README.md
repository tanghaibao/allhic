# ALLHIC: Genome scaffolding based on HiC data

[![Travis-CI](https://travis-ci.org/tanghaibao/allhic.svg?branch=master)](https://travis-ci.org/tanghaibao/allhic)

**This software is currently under active development. DO NOT USE.**

## Installation

The easiest way to install allhic is to download the latest binary from
the [releases](https://github.com/tanghaibao/allhic/releases) and make sure to
`chmod +x` the resulting binary.

If you are using [go](https://github.com/golang/go), you can build from source with:
```
go get -u -t -v github.com/tanghaibao/allhic/...
go install github.com/tanghaibao/allhic/cmd/allhic
```

## Usage

### Prune

Prune bamfile to remove weak links. WIP.

### Partition

Given a target `k`, number of partitions, the goal of the partitioning
is to separate all the contigs into separate clusters. As with all
clustering algorithm, there is an optimization goal here. The
LACHESIS algorithm is a hierarchical clustering algorithm using
average links. ALLHIC uses a community detection method based on
[Newman 2006](http://www.pnas.org/content/103/23/8577.full),
using eigen decomposition of the modularity matrix.

![networkbefore](script/graph.png)
![networkafter](script/graph.partitioned.png)

```bash
allhic partition tests/test.bam
```

### Optimize

Given a set of Hi-C contacts between contigs, as specified in the
clmfile, reconstruct the highest scoring ordering and orientations
for these contigs.

Optimize uses Genetic Algorithm (GA) to search for the best scoring solution.

![ga](tests/test-movie.gif)

```bash
allhic optimize tests/test.clm
```

### Build

Build genome release. WIP.
