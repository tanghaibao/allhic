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
	Name        string
	Clmfile     string
	Idsfile     string
	tigToSize   map[string]int
	tigActive   map[string]bool
	tigToIdx    map[string]int
	activeSizes []int
	contacts    []Contact
}

// Contact stores how many links between two contigs
type Contact struct {
	a      string
	b      string
	ao     int8
	bo     int8
	nlinks int
	dists  []int
}

// InitCLMFile is the constructor for CLMFile
func InitCLMFile(Clmfile string) *CLMFile {
	p := new(CLMFile)
	p.Name = RemoveExt(path.Base(Clmfile))
	p.Clmfile = Clmfile
	p.Idsfile = RemoveExt(Clmfile) + ".ids"
	p.tigToSize = make(map[string]int)
	p.tigActive = make(map[string]bool)
	p.tigToIdx = make(map[string]int)

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
		r.tigActive[tig] = true
	}

	fmt.Println(r.tigToSize)
}

// Map orientations to ints
func ff(c byte) int8 {
	if c == '-' {
		return -1
	}
	return 1
}
func rr(c byte) int8 {
	if c == '-' {
		return 1
	}
	return -1
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
		var contact = Contact{at, bt, ff(ao), ff(bo), nlinks, dists}
		r.contacts = append(r.contacts, contact)
		// strandedness := ao == bo
		// fmt.Println(contact.a, contact.b, contact.nlinks, strandedness)
	}
}

// UpdateTigToIdx maps tigs to indices in the current active tigs
func (r *CLMFile) UpdateTigToIdx() {
	idx := 0
	r.activeSizes = make([]int, len(r.tigActive))
	for k := range r.tigActive {
		r.tigToIdx[k] = idx
		r.activeSizes[idx] = r.tigToSize[k]
		idx++
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
func (r *CLMFile) Activate() []int {
	tour := make([]int, len(r.tigActive))
	r.UpdateTigToIdx()
	r.reportActive()
	for i := 0; i < len(r.tigActive); i++ {
		tour[i] = i
	}

	return tour
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
	N := len(r.tigActive)
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
