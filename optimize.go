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
	"math/rand"
	"os"
	"strings"
)

// Optimizer runs the order-and-orientation procedure, given a clmfile
type Optimizer struct {
	REfile    string
	Clmfile   string
	RunGA     bool
	Resume    bool
	Seed      int64
	NPop      int
	NGen      int
	MutProb   float64
	CrossProb float64
	rng       *rand.Rand
	// Output files
	OutTourFile string
}

// Run kicks off the Optimizer
func (r *Optimizer) Run() error {
	r.rng = rand.New(rand.NewSource(r.Seed))
	clm, err := NewCLM(r.Clmfile, r.REfile)
	if err != nil {
		return err
	}
	tourfile := RemoveExt(r.REfile) + ".tour"

	// Load tourfile if it exists
	if _, err := os.Stat(tourfile); r.Resume && err == nil {
		log.Noticef("Found existing tour file `%s`", tourfile)
		err = clm.parseTourFile(tourfile)
		if err != nil {
			return err
		}
		// Rename the tour file
		backupTourFile := tourfile + ".sav"
		err = os.Rename(tourfile, backupTourFile)
		if err != nil {
			return err
		}
		log.Noticef("Backup `%s` to `%s`", tourfile, backupTourFile)
	}

	clm.Activate(true, r.rng)

	// tourfile logs the intermediate configurations
	log.Noticef("Optimization history logged to `%s`", tourfile)
	fwtour, _ := os.Create(tourfile)
	r.OutTourFile = tourfile

	err = clm.printTour(os.Stdout, clm.Tour, "INIT")
	if err != nil {
		return err
	}
	err = clm.printTour(fwtour, clm.Tour, "INIT")
	if err != nil {
		return err
	}

	if r.RunGA {
		for phase := 1; phase < 3; phase++ {
			err = clm.OptimizeOrdering(fwtour, r, phase)
			if err != nil {
				return err
			}
		}
	}

	for phase := 1; ; phase++ {
		tag1, tag2, err := clm.OptimizeOrientations(fwtour, phase)
		if err != nil {
			return err
		}
		if tag1 == REJECT && tag2 == REJECT {
			log.Noticef("Terminating ... no more %v", ACCEPT)
			break
		}
	}
	err = clm.printTour(os.Stdout, clm.Tour, "FINAL")
	if err != nil {
		return err
	}
	log.Notice("Success")
	return fwtour.Close()
}

// OptimizeOrdering changes the ordering of contigs by Genetic Algorithm
func (r *CLM) OptimizeOrdering(fwtour *os.File, opt *Optimizer, phase int) error {
	_, err := r.GARun(fwtour, opt, phase)
	return err
}

// OptimizeOrientations changes the orientations of contigs by using heuristic flipping algorithms.
func (r *CLM) OptimizeOrientations(fwtour *os.File, phase int) (string, string, error) {
	tag1 := r.flipWhole()
	err := r.printTour(fwtour, r.Tour, fmt.Sprintf("FLIPWHOLE%d", phase))
	if err != nil {
		return "", "", err
	}
	tag2 := r.flipOne()
	err = r.printTour(fwtour, r.Tour, fmt.Sprintf("FLIPONE%d", phase))
	return tag1, tag2, err
}

// parseTourFile parses tour file
// Only the last line is retained anc onverted into a Tour
func parseTourFile(filename string) ([]string, error) {
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	log.Noticef("Parse tour file `%s`", filename)

	reader := bufio.NewReader(f)
	var words []string
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		if row[0] == '>' { // header
			continue
		}
		words = strings.Split(row, " ")
	}
	err = f.Close()
	return words, err
}

// prepareTour prepares a boilerplate for an empty tour
func (r *CLM) prepareTour() {
	r.Signs = make([]byte, len(r.Tigs))
	for _, tig := range r.Tigs {
		tig.IsActive = false
	}
}

// parseTourFile parses tour file
// Only the last line is retained anc onverted into a Tour
func (r *CLM) parseTourFile(filename string) error {
	words, err := parseTourFile(filename)
	if err != nil {
		return err
	}
	r.prepareTour()

	tigs := make([]Tig, 0)
	for _, word := range words {
		tigName, tigOrientation := word[:len(word)-1], word[len(word)-1]
		idx, ok := r.tigToIdx[tigName]
		if !ok {
			log.Errorf("Contig %s not found!", tigName)
			continue
		}
		tigs = append(tigs, Tig{
			Idx: idx,
		})
		r.Signs[idx] = tigOrientation
		r.Tigs[idx].IsActive = true
	}
	r.Tour.Tigs = tigs
	return r.printTour(os.Stdout, r.Tour, "INIT")
}

// parseClustersFile parses clusters file
func (r *CLM) parseClustersFile(clustersfile string, group int) error {
	recs, err := ReadCSVLines(clustersfile)
	if err != nil {
		return err
	}
	r.prepareTour()

	rec := recs[group]
	names := strings.Split(rec[2], " ")
	tigs := make([]Tig, 0)
	for _, tigName := range names {
		idx, ok := r.tigToIdx[tigName]
		if !ok {
			log.Errorf("Contig %s not found!", tigName)
			continue
		}
		tigs = append(tigs, Tig{
			Idx: idx,
		})
		r.Signs[idx] = '+'
		r.Tigs[idx].IsActive = true
	}
	r.Tour.Tigs = tigs
	return r.printTour(os.Stdout, r.Tour, "INIT")
}

// printTour logs the current tour to file
func (r *CLM) printTour(fwtour *os.File, tour Tour, label string) error {
	_, err := fwtour.WriteString(">" + label + "\n")
	if err != nil {
		return err
	}
	atoms := make([]string, tour.Len())
	for i := 0; i < tour.Len(); i++ {
		idx := tour.Tigs[i].Idx
		atoms[i] = r.Tigs[idx].Name + string(r.Signs[idx])
	}
	_, err = fwtour.WriteString(strings.Join(atoms, " ") + "\n")
	return err
}
