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
	"os"
	"strings"
)

// Optimizer runs the order-and-orientation procedure, given a clmfile
type Optimizer struct {
	Clmfile string
}

// Run kicks off the Optimizer
func (r *Optimizer) Run() {
	clm := InitCLMFile(r.Clmfile)
	tourfile := RemoveExt(r.Clmfile) + ".tour"
	fwtour, _ := os.Create(tourfile)
	defer fwtour.Close()

	shuffle := true
	clm.Activate(shuffle)
	clm.OptimizeOrdering(fwtour)
	clm.OptimizeOrientations(fwtour)
}

// OptimizeOrdering changes the ordering of contigs by Genetic Algorithm
func (r *CLMFile) OptimizeOrdering(fwtour *os.File) {
	gaTour := GARun(r.Tour, 100, 2000, .2)
	r.Tour = gaTour
	r.pruneTour()
	r.PrintTour(fwtour, r.Tour, "test")
}

// OptimizeOrientations changes the orientations of contigs by using heuristic flipping algorithms.
func (r *CLMFile) OptimizeOrientations(fwtour *os.File) {
	r.flipAll()
}

// PrintTour logs the current tour to file
func (r *CLMFile) PrintTour(fwtour *os.File, tour Tour, label string) {
	fwtour.WriteString(">" + label + "\n")
	atoms := make([]string, tour.Len())
	for i := 0; i < tour.Len(); i++ {
		atoms[i] = r.Tigs[tour.Tigs[i].Idx].Name
	}
	fwtour.WriteString(strings.Join(atoms, " ") + "\n")
}
