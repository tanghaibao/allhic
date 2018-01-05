/**
 * Filename: /Users/bao/code/allhic/allhic/clm.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Monday, January 1st 2018, 5:57:00 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"bufio"
	"fmt"
	"os"
	"path"
	"strconv"
	"strings"
)

// CLMFile has the following format:
//
// tig00046211+ tig00063795+       1       53173
// tig00046211+ tig00063795-       1       116050
// tig00046211- tig00063795+       1       71155
// tig00046211- tig00063795-       1       134032
// tig00030676+ tig00077819+       7       136407 87625 87625 106905 102218 169660 169660
// tig00030676+ tig00077819-       7       126178 152952 152952 35680 118923 98367 98367
// tig00030676- tig00077819+       7       118651 91877 91877 209149 125906 146462 146462
// tig00030676- tig00077819-       7       108422 157204 157204 137924 142611 75169 75169
type CLMFile struct {
	Name                 string
	Clmfile              string
	Idsfile              string
	tigToSize            map[string]int
	tigToIdx             map[string]int
	activeTigs           []string
	activeSizes          []int
	contacts             []Contact                // Array of contacts
	orientations         map[Pair]OrientedContact // (tigA, tigB) => strandedness x nlinks
	contactsOrientations map[Pair][4]GArray       // (tigA, tigB) => [0..3]gdists, 0..3 are orientations
}

// Pair contains two contigs in contact
type Pair struct {
	a string
	b string
}

// Contact stores how many links between two contigs
type Contact struct {
	a      string
	b      string
	nlinks int
	dists  []int
}

// OrientedContact stores only one configuration per pair of tigs
type OrientedContact struct {
	strandedness int
	nlinks       int
	meanDist     int
}

// Tig stores the index to activeTigs and size of the tig
type Tig struct {
	Idx  int
	Size int
}

// Tour stores a number of tigs along with 2D matrices for evaluation
type Tour struct {
	Tigs []Tig
	M    [][]int
}

// InitCLMFile is the constructor for CLMFile
func InitCLMFile(Clmfile string) *CLMFile {
	p := new(CLMFile)
	p.Name = RemoveExt(path.Base(Clmfile))
	p.Clmfile = Clmfile
	p.Idsfile = RemoveExt(Clmfile) + ".ids"
	p.tigToSize = make(map[string]int)
	p.tigToIdx = make(map[string]int)
	p.orientations = make(map[Pair]OrientedContact)
	p.contactsOrientations = make(map[Pair][4]GArray)

	p.ParseIds()
	p.ParseClm()

	return p
}

// ParseIds parses the idsfile into data stored in CLMFile.
// IDS file has a list of contigs that need to be ordered. 'recover',
// keyword, if available in the third column, is less confident.
// tig00015093     46912
// tig00035238     46779   recover
// tig00030900     119291
func (r *CLMFile) ParseIds() {
	file, _ := os.Open(r.Idsfile)
	log.Noticef("Parse idsfile `%s`", r.Idsfile)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		tig := words[0]
		size, _ := strconv.Atoi(words[1])
		r.tigToSize[tig] = size
		r.activeTigs = append(r.activeTigs, tig)
	}

	fmt.Println(r.tigToSize)
}

// ff map orientations to bit ('+' => 0, '-' => 1)
func ff(c byte) (bit byte) {
	if c != '-' {
		bit = 1
	}
	return
}

// bb map two orientations to a number 0..3
func bb(a, b byte) byte {
	return a<<1 + b
}

// ParseClm parses the clmfile into data stored in CLMFile.
func (r *CLMFile) ParseClm() {
	file, _ := os.Open(r.Clmfile)
	log.Noticef("Parse clmfile `%s`", r.Clmfile)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		row := strings.TrimSpace(scanner.Text())
		words := strings.Split(row, "\t")
		abtig := strings.Split(words[0], " ")
		atig, btig := abtig[0], abtig[1]
		at, ao := atig[:len(atig)-1], atig[len(atig)-1]
		bt, bo := btig[:len(btig)-1], btig[len(btig)-1]

		// Make sure both contigs are in the ids file
		if _, ok := r.tigToSize[at]; !ok {
			continue
		}
		if _, ok := r.tigToSize[bt]; !ok {
			continue
		}

		nlinks, _ := strconv.Atoi(words[1])
		// Convert all distances to int
		var dists []int
		for _, dist := range strings.Split(words[2], " ") {
			d, _ := strconv.Atoi(dist)
			dists = append(dists, d)
		}

		// Store all these info in contacts
		contact := Contact{at, bt, nlinks, dists}
		// gdists := GoldenArray(dists)
		meanDist := HmeanInt(dists, GRLB, GRUB)
		// fmt.Println(at, bt, dists, gdists)
		strandedness := 1
		if ao != bo {
			strandedness = -1
		}
		r.contacts = append(r.contacts, contact)
		pair := Pair{at, bt}
		r.orientations[pair] = OrientedContact{strandedness, nlinks, meanDist}
		// r.contactsOrientations[pair][bb(ff(ao), ff(bo))] = gdists
		// r.contactsOrientations[pair][bb()]
		// strandedness := ao == bo
	}
}

// UpdateTigToIdx maps tigs to indices in the current active tigs
func (r *CLMFile) UpdateTigToIdx() {
	r.activeSizes = make([]int, len(r.activeTigs))
	for i, k := range r.activeTigs {
		r.tigToIdx[k] = i
		r.activeSizes[i] = r.tigToSize[k]
	}
}

// Activate selects contigs in the current partition. This is the setup phase of the
// algorithm, and supports two modes:
// - "de novo": This is useful at the start of a new run where no tours are
//    available. We select the strong contigs that have significant number
//    of links to other contigs in the partition. We build a histogram of
//    link density (# links per bp) and remove the contigs that appear to be
//    outliers. The orientations are derived from the matrix decomposition
//    of the pairwise strandedness matrix O.
// - "hotstart": This is useful when there was a past run, with a given
//    tourfile. In this case, the active contig list and orientations are
//    derived from the last tour in the file.
func (r *CLMFile) Activate() (tour Tour) {
	tour.Tigs = make([]Tig, len(r.activeTigs))
	r.UpdateTigToIdx()
	r.reportActive()
	for i := 0; i < len(r.activeTigs); i++ {
		tour.Tigs[i] = Tig{i, r.activeSizes[i]}
	}
	tour.M = r.M()

	return
}

// reportActive prints out a quick message on number of active tigs
func (r *CLMFile) reportActive() {
	sumLength := 0
	for _, v := range r.activeSizes {
		sumLength += v
	}
	log.Noticef("Active contigs: %d (length=%d)",
		len(r.activeSizes), sumLength)
}

// M yields a contact frequency matrix, where each cell contains how many
// links between i-th and j-th contig
func (r *CLMFile) M() [][]int {
	N := len(r.activeTigs)
	P := make([][]int, N)
	for i := 0; i < N; i++ {
		P[i] = make([]int, N)
	}

	for _, contact := range r.contacts {
		ai, aok := r.tigToIdx[contact.a]
		if !aok {
			continue
		}
		bi, bok := r.tigToIdx[contact.b]
		if !bok {
			continue
		}
		P[ai][bi] = contact.nlinks
		P[bi][ai] = contact.nlinks
	}

	return P
}
