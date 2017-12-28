package main

import (
	"fmt"

	"./allhic"

	"github.com/docopt/docopt-go"
)

func main() {
	usage := `ALLHIC: genome scaffolding based on Hi-C data

Usage:
  allhic partition
  allhic optimize

Options:
  -h --help       Show this screen.
  --version       Show version.`

	arguments, _ := docopt.Parse(usage, nil, true, "ALLHIC 0.7.11", false)
	fmt.Println(arguments)

	p := allhic.Partitioner{"tests/prunning.sub.bam"}
	fmt.Println(p.Bamfile)
	p.CountLinks()
}
