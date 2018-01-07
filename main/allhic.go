/**
 * Filename: /Users/htang/code/allhic/main.go
 * Path: /Users/htang/code/allhic
 * Created Date: Wednesday, January 3rd 2018, 11:21:45 am
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

package main

import (
	"fmt"
	"os"

	"github.com/tanghaibao/allhic"
	"github.com/docopt/docopt-go"
	logging "github.com/op/go-logging"
)

// main is the entrypoint for the entire program, routes to commands
func main() {
	usage := `ALLHIC: genome scaffolding based on Hi-C data

Usage:
  allhic partition [options] <bamfile>
  allhic optimize [options] <clmfile>`

	logging.SetBackend(allhic.BackendFormatter)
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, usage)
		return
	}

	command := os.Args[1]
	if command == "partition" {
		partitionMain()
	} else if command == "optimize" {
		optimizeMain()
	} else {
		fmt.Fprintln(os.Stderr, usage)
		return
	}
}

// partitionMain is the entrypoint for partition
func partitionMain() {
	usage := `ALLHIC: genome scaffolding based on Hi-C data

Partition function:
Given a target k, number of partitions, the goal of the partitioning is to
separate all the contigs into separate clusters. As with all clustering
algorithm, there is an optimization goal here. The LACHESIS algorithm is
a hierarchical clustering algorithm using average links.

Usage:
    allhic partition <bamfile> [options]

Options:
    --help       Show this screen.
    --version    Show version.`

	args, _ := docopt.Parse(usage, nil, true, allhic.Version, false)
	p := allhic.Partitioner{args["<bamfile>"].(string)}
	p.Run()
}

// optimizeMain is the entrypoint for optimize
func optimizeMain() {
	usage := `ALLHIC: genome scaffolding based on Hi-C data

Optimize function:
Given a set of Hi-C contacts between contigs, as specified in the
clmfile, reconstruct the highest scoring ordering and orientations
for these contigs.

Usage:
    allhic optimize <clmfile> [options]

Options:
	--help       Show this screen
	--version    Show version
	--skipGA     Skip GA step`

	args, _ := docopt.Parse(usage, nil, true, allhic.Version, false)
	runGA := !(args["--skipGA"].(bool))
	p := allhic.Optimizer{args["<clmfile>"].(string), runGA}
	p.Run()
}
