# ALLHIC: Genome scaffolding based on HiC data

         _       _____     _____     ____  ____  _____   ______
        / \     |_   _|   |_   _|   |_   ||   _||_   _|.' ___  |
       / _ \      | |       | |       | |__| |    | | / .'   \_|
      / ___ \     | |   _   | |   _   |  __  |    | | | |
    _/ /   \ \_  _| |__/ | _| |__/ | _| |  | |_  _| |_\ `.___.'\
    |____| |____||________||________||____||____||_____|`.____ .'

[![Travis-CI](https://travis-ci.org/tanghaibao/allhic.svg?branch=master)](https://travis-ci.org/tanghaibao/allhic)
[![GOreport](https://goreportcard.com/badge/github.com/tanghaibao/allhic)](https://goreportcard.com/report/github.com/tanghaibao/allhic)

**This software is currently under active development. DO NOT USE.**

## Installation

The easiest way to install allhic is to download the latest binary from
the [releases](https://github.com/tanghaibao/allhic/releases) and make sure to
`chmod +x` the resulting binary.

If you are using [go](https://github.com/golang/go), you can build from source with:

```console
go get -u -t -v github.com/tanghaibao/allhic/...
go install github.com/tanghaibao/allhic/cmd/allhic
```

## Usage

### Prune

Prune bamfile to remove weak links. WIP.

### Extract

Extract does a fair amount of preprocessing: 1) extract inter-contig links into a more compact form, specifically into `.clm`; 2) extract intra-contig links and build a distribution; 3) count up the restriction sites to be used in normalization (similar to LACHESIS); 4) bundles the inter-contig links into pairs of contigs.

```console
allhic extract tests/test.bam tests/test.fasta
```

### Partition

Given a target `k`, number of partitions, the goal of the partitioning
is to separate all the contigs into separate clusters. As with all
clustering algorithm, there is an optimization goal here. The
LACHESIS algorithm is a hierarchical clustering algorithm using
average links. ALLHIC uses a community detection method based on
[Newman 2006](http://www.pnas.org/content/103/23/8577.full),
using eigen decomposition of the modularity matrix.

![networkbefore](script/graph-s.png)
![networkafter](script/graph-s.partitioned.png)

```console
allhic partition tests/test.counts_GATC.txt tests/test.pairs.txt
```

### Optimize

Given a set of Hi-C contacts between contigs, as specified in the
clmfile, reconstruct the highest scoring ordering and orientations
for these contigs.

Optimize uses Genetic Algorithm (GA) to search for the best scoring solution.

![ga](tests/test-movie.gif)

```console
allhic optimize tests/test.counts_GATC.txt tests/test.clm
```

### Build

Build genome release. WIP.

## Pipeline

Following the 4 steps of `prune`, `extract`, `partition`, `optimize`

```console
allhic extract T4_Chr1/{prunning.sub.bam,seq.fasta}
allhic partition T4_Chr1/{prunning.sub.counts_GATC.txt,prunning.sub.pairs.txt} 2
allhic optimize T4_Chr1/{prunning.sub.counts_GATC.txt,prunning.sub.clm}
allhic build T4_Chr/{prunning.sub.tour,seq.fasta}
```

## WIP features

- [x] Add restriction enzyme for better normalization of contig lengths
- [ ] Add test suites
