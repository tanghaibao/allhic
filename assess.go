/*
 * Filename: /Users/htang/code/allhic/assess.go
 * Path: /Users/htang/code/allhic
 * Created Date: Tuesday, June 19th 2018, 4:34:11 pm
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"os"

	"github.com/biogo/hts/bam"
)

// Assesser takes input of bamfile and bedfile and output per contig confidence
// in the orientation
//
// Summary of algorithm:
// Step 1. Take all intra-contig links and build the background distribution
// Step 2. Loop through each contig, compute the likelihood of all links coming
//         out of the contig, assuming + orientation, and - orientation, separately
// Step 3. Normalize the likelihood to get the posterior probability (implicit assumption)
//         of equal prior probability for each contig
type Assesser struct {
	Bamfile string
}

// Run calls the Assessor
func (r *Assesser) Run() {
	r.ExtractIntraContigLinks()
	r.ComputeLikelihood()
	r.ComputePosteriorProb()
}

// ExtractIntraContigLinks builds the probability distribution of link sizes
func (r *Assesser) ExtractIntraContigLinks() {
	fh, _ := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()
}

// ComputeLikelihood computes the likelihood of link sizes assuming + orientation
// and - orientation, respectively
func (r *Assesser) ComputeLikelihood() {

}

// ComputePosteriorProb computes the posterior probability of the orientations
func (r *Assesser) ComputePosteriorProb() {

}
