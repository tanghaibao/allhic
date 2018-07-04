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
	"bufio"
	"fmt"
	"io"
	"os"
	"strings"
)

// Optimizer runs the order-and-orientation procedure, given a clmfile
type Optimizer struct {
	Clmfile   string
	RunGA     bool
	StartOver bool
	Seed      int64
	NPop      int
	NGen      int
	MutProb   float64
	CrossProb float64
}

// Run kicks off the Optimizer
func (r *Optimizer) Run() {
	clm := NewCLM(r.Clmfile)
	tourfile := RemoveExt(r.Clmfile) + ".tour"
	var shuffle bool

	// Load tourfile if it exists
	if _, err := os.Stat(tourfile); !r.StartOver && err == nil {
		log.Noticef("Found existing tour file")
		clm.parseTourFile(tourfile)
		clm.PrintTour(os.Stdout, clm.Tour, "INIT")
		// Rename the tour file
		backupTourFile := tourfile + ".sav"
		os.Rename(tourfile, backupTourFile)
		log.Noticef("Backup `%s` to `%s`", tourfile, backupTourFile)
	} else {
		shuffle = true
	}

	clm.Activate(shuffle)

	// tourfile logs the intermediate configurations
	fwtour, _ := os.Create(tourfile)
	defer fwtour.Close()
	clm.PrintTour(fwtour, clm.Tour, "INIT")

	if r.RunGA {
		for phase := 1; phase < 3; phase++ {
			clm.OptimizeOrdering(fwtour, r, phase)
		}
	}

	for phase := 1; ; phase++ {
		tag1, tag2 := clm.OptimizeOrientations(fwtour, phase)
		if tag1 == REJECT && tag2 == REJECT {
			log.Noticef("Terminating ... no more %v", ACCEPT)
			break
		}
	}
	log.Notice("Success")
}

// OptimizeOrdering changes the ordering of contigs by Genetic Algorithm
func (r *CLM) OptimizeOrdering(fwtour *os.File, opt *Optimizer, phase int) {
	r.GARun(fwtour, opt, phase)
	r.pruneTour()
}

// OptimizeOrientations changes the orientations of contigs by using heuristic flipping algorithms.
func (r *CLM) OptimizeOrientations(fwtour *os.File, phase int) (string, string) {
	tag1 := r.flipWhole()
	r.PrintTour(fwtour, r.Tour, fmt.Sprintf("FLIPWHOLE%d", phase))
	tag2 := r.flipOne()
	r.PrintTour(fwtour, r.Tour, fmt.Sprintf("FLIPONE%d", phase))
	return tag1, tag2
}

// parseTourFile parses tour file
// Only the last line is retained anc onverted into a Tour
func (r *CLM) parseTourFile(filename string) {
	f, _ := os.Open(filename)
	log.Noticef("Parse tour file `%s`", filename)
	defer f.Close()

	reader := bufio.NewReader(f)
	var tigs []Tig
	r.Signs = make([]byte, len(r.Tigs))
	for _, tig := range r.Tigs {
		tig.IsActive = false
	}
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		if row[0] == '>' { // header
			continue
		}
		words := strings.Split(row, " ")
		tigs = make([]Tig, len(words))
		for i, word := range words {
			tigName, tigOrientation := word[:len(word)-1], word[len(word)-1]
			idx, ok := r.tigToIdx[tigName]
			if !ok {
				log.Errorf("Contig %s not found!", tigName)
			}
			tigs[i].Idx = idx
			r.Signs[idx] = tigOrientation
			r.Tour.Tigs = tigs
			r.Tigs[idx].IsActive = true
		}
	}
}

// PrintTour logs the current tour to file
func (r *CLM) PrintTour(fwtour *os.File, tour Tour, label string) {
	fwtour.WriteString(">" + label + "\n")
	atoms := make([]string, tour.Len())
	for i := 0; i < tour.Len(); i++ {
		idx := tour.Tigs[i].Idx
		atoms[i] = r.Tigs[idx].Name + string(r.Signs[idx])
	}
	fwtour.WriteString(strings.Join(atoms, " ") + "\n")
}
