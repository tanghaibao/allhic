package main

import (
	"./allhic"

	"github.com/docopt/docopt-go"
	logging "github.com/op/go-logging"
)

func main() {
	usage := `ALLHIC: genome scaffolding based on Hi-C data

Usage:
  allhic partition
  allhic optimize

Options:
  -h --help       Show this screen.
  --version       Show version.`

	arguments, _ := docopt.Parse(usage, nil, true, "ALLHIC 0.8.1", false)
	//fmt.Println(arguments)

	logging.SetBackend(allhic.BackendFormatter)

	if arguments["partition"].(bool) {
		p := allhic.Partitioner{"tests/prunning.sub.bam"}
		p.CountLinks()
	} else if arguments["optimize"].(bool) {
		p := allhic.Optimizer{"tests/test.clm"}
		p.Run()
	}
}
