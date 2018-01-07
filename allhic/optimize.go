/**
 * Filename: /Users/bao/code/allhic/allhic/optimize.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Tuesday, January 2nd 2018, 10:00:33 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"fmt"
	"os"
	"strings"
)

// Optimizer runs the order-and-orientation procedure, given a clmfile
type Optimizer struct {
	Clmfile string
	RunGA   bool
}

// Run kicks off the Optimizer
func (r *Optimizer) Run() {
	clm := InitCLMFile(r.Clmfile)
	tourfile := RemoveExt(r.Clmfile) + ".tour"

	shuffle := true
	clm.Activate(shuffle)

	// tourfile logs the intermediate configurations
	fwtour, _ := os.Create(tourfile)
	defer fwtour.Close()
	clm.PrintTour(fwtour, clm.Tour, "INIT")

	if r.RunGA {
		for phase := 1; phase < 3; phase++ {
			clm.OptimizeOrdering(fwtour, phase)
		}
	}

	for phase := 1; ; phase++ {
		tag1, tag2 := clm.OptimizeOrientations(fwtour, phase)
		if tag1 == REJECT && tag2 == REJECT {
			log.Noticef("Terminating ... no more %v", ACCEPT)
			break
		}
	}
}

// OptimizeOrdering changes the ordering of contigs by Genetic Algorithm
func (r *CLMFile) OptimizeOrdering(fwtour *os.File, phase int) {
	gaTour := GARun(r.Tour, 100, 2000, .2)
	r.Tour = gaTour
	r.pruneTour()
	r.PrintTour(fwtour, r.Tour, fmt.Sprintf("GA%d", phase))
}

// OptimizeOrientations changes the orientations of contigs by using heuristic flipping algorithms.
func (r *CLMFile) OptimizeOrientations(fwtour *os.File, phase int) (string, string) {
	tag1 := r.flipWhole()
	r.PrintTour(fwtour, r.Tour, fmt.Sprintf("FLIPWHOLE%d", phase))
	tag2 := r.flipOne()
	r.PrintTour(fwtour, r.Tour, fmt.Sprintf("FLIPONE%d", phase))
	return tag1, tag2
}

// PrintTour logs the current tour to file
func (r *CLMFile) PrintTour(fwtour *os.File, tour Tour, label string) {
	fwtour.WriteString(">" + label + "\n")
	atoms := make([]string, tour.Len())
	for i := 0; i < tour.Len(); i++ {
		idx := tour.Tigs[i].Idx
		atoms[i] = r.Tigs[idx].Name + string(r.Signs[idx])
	}
	fwtour.WriteString(strings.Join(atoms, " ") + "\n")
}
