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
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

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
	Bamfile             string
	Bedfile             string
	Seqid               string
	contigs             []BedLine
	contigIntraLinks    []int   // Contig link sizes for all intra-contig links
	contigInterLinksFwd [][]int // Contig link sizes assuming same dir
	contigInterLinksRev [][]int // Contig link sizes assuming other dir
}

// BedLine stores the information from each line in the bedfile
type BedLine struct {
	seqid string
	start int
	end   int
	name  string
}

// Run calls the Assessor
func (r *Assesser) Run() {
	r.ReadBed()
	r.ExtractContigLinks()
	r.ComputeLikelihood()
	r.ComputePosteriorProb()
}

// ReadBed parses the bedfile to extract the start and stop for all the contigs
func (r *Assesser) ReadBed() {
	fh, _ := os.Open(r.Bedfile)
	log.Noticef("Parse bedfile `%s`", r.Bedfile)
	reader := bufio.NewReader(fh)

	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		words := strings.Split(row, "\t")
		seqid := words[0]
		if seqid != r.Seqid {
			continue
		}
		start, _ := strconv.Atoi(words[1])
		end, _ := strconv.Atoi(words[2])
		r.contigs = append(r.contigs, BedLine{
			seqid: seqid,
			start: start - 1,
			end:   end,
			name:  words[3],
		})
	}

	sort.Slice(r.contigs, func(i, j int) bool {
		return r.contigs[i].start < r.contigs[j].start
	})
	log.Noticef("A total of %d contigs imported", len(r.contigs))
}

// checkInRange checks if a point position is within range
func checkInRange(pos, start, end int) bool {
	return start <= pos && pos < end
}

// ExtractContigLinks builds the probability distribution of link sizes
func (r *Assesser) ExtractContigLinks() {
	fh, _ := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	// Import links into pairs of contigs
	r.contigInterLinksFwd = make([][]int, len(r.contigs))
	r.contigInterLinksRev = make([][]int, len(r.contigs))
	var a, b int
	nIntraLinks := 0
	nInterLinks := 0
	ci := 0 // Use this to index into r.contigs, the current contig under consideration
	for {
		rec, err := br.Read()
		if err != nil {
			if err != io.EOF {
				log.Error(err)
			}
			break
		}

		// Restrict the links to be within the current chromosome
		at, bt := rec.Ref.Name(), rec.MateRef.Name()
		if at != r.Seqid || bt != r.Seqid {
			continue
		}

		//         read1                                               read2
		//     ---a-- X|----- dist = a2 ----|         |--- dist = b ---|X ------ b2 ------
		//     ==============================         ====================================
		//             C1 (length L1)       |----D----|         C2 (length L2)
		a, b = rec.Pos, rec.MatePos
		if a < r.contigs[ci].start {
			continue
		}

		// Now we need to check if this pair of positions is a intra-contig or inter-contig link
		// If the intervals are disjoint and the mapping lies between the intervals, then this could
		// lead to a problem
		for a > r.contigs[ci].end {
			ci++
			fmt.Println(r.contigs[ci], a, nIntraLinks, nInterLinks)
		}

		link := abs(a - b)
		// For intra-contig link it's easy, just store the distance between two ends
		// An intra-contig link
		if checkInRange(b, r.contigs[ci].start, r.contigs[ci].end) {
			r.contigIntraLinks = append(r.contigIntraLinks, link)
			nIntraLinks++
			continue
		}

		// For inter-contig link this is a bit tricky:
		// - for the same direction as in the AGP/BED, the distance is the real distance
		// - for the opposite direction as in the AGP/BED, we need to flip the contig,
		//   which means, contig_start + contig_end - link_start
		//   To check this is correct, link_start = contig_start ==> contig_end
		//                        and, link_start = contig_end => contig_start
		// An inter-contig link
		r.contigInterLinksFwd[ci] = append(r.contigInterLinksFwd[ci], link)
		// Assuming flipped orientation
		link = abs(r.contigs[ci].start + r.contigs[ci].end - a - b)
		r.contigInterLinksRev[ci] = append(r.contigInterLinksRev[ci], link)
		nInterLinks++
	}
	log.Noticef("A total of %d intra-contig links and %d inter-contig links imported",
		nIntraLinks, nInterLinks)
}

// ComputeLikelihood computes the likelihood of link sizes assuming + orientation
// and - orientation, respectively
func (r *Assesser) ComputeLikelihood() {

}

// ComputePosteriorProb computes the posterior probability of the orientations
func (r *Assesser) ComputePosteriorProb() {

}
