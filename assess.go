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
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
)

// probCutoff is the minimum level of prob required
const probCutoff = .95

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
	Bamfile       string
	Bedfile       string
	Seqid         string
	seq           *ContigInfo
	model         *LinkDensityModel
	contigs       []BedLine
	interLinksFwd [][]int // Contig link sizes assuming same dir
	interLinksRev [][]int // Contig link sizes assuming other dir
	postprob      []float64
}

// BedLine stores the information from each line in the bedfile
type BedLine struct {
	seqid string
	start int
	end   int
	name  string
	size  int
}

// Run calls the Assessor
func (r *Assesser) Run() {
	r.readBed()
	r.extractContigLinks()
	r.makeModel(r.Seqid + ".distribution.txt")
	r.computePosteriorProb()
	r.writePostProb(r.Seqid + ".postprob.txt")
	log.Notice("Success")
}

// makeModel computes the norms and bins separately to derive an empirical link size
// distribution, then power law is inferred for extrapolating higher values
func (r *Assesser) makeModel(outfile string) {
	contigSizes := []int{}
	for _, contig := range r.contigs {
		contigSizes = append(contigSizes, contig.size)
	}
	m := NewLinkDensityModel()
	m.makeBins()
	m.makeNorms(contigSizes)
	m.countBinDensities([]*ContigInfo{r.seq})
	m.writeDistribution(outfile)
	r.model = m
}

// writePostProb writes the final posterior probability to file
func (r *Assesser) writePostProb(outfile string) {
	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()

	fmt.Fprintf(w, "#SeqID\tStart\tEnd\tContig\tPostProb\n")
	for i, contig := range r.contigs {
		fmt.Fprintf(w, "%s\t%d\t%d\t%s\t%.4f\n",
			contig.seqid, contig.start, contig.end, contig.name, r.postprob[i])
	}

	w.Flush()
	log.Noticef("Posterior probability written to `%s`", outfile)
}

// readBed parses the bedfile to extract the start and stop for all the contigs
func (r *Assesser) readBed() {
	fh := mustOpen(r.Bedfile)
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
		// start-- // To handle sometimes 1-based offset
		end, _ := strconv.Atoi(words[2])
		r.contigs = append(r.contigs, BedLine{
			seqid: seqid,
			start: start,
			end:   end,
			name:  words[3],
			size:  end - start,
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

// extractContigLinks builds the probability distribution of link sizes
func (r *Assesser) extractContigLinks() {
	fh := mustOpen(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	// We need the size of the SeqId to compute expected number of links
	var s *ContigInfo
	refs := br.Header().Refs()
	for _, ref := range refs {
		if ref.Name() == r.Seqid {
			s = &ContigInfo{
				name:   ref.Name(),
				length: ref.Len(),
				links:  []int{},
			}
			r.seq = s
			break
		}
	}

	log.Noticef("Seq `%s` has size %d", s.name, s.length)

	// Import links into pairs of contigs
	r.interLinksFwd = make([][]int, len(r.contigs))
	r.interLinksRev = make([][]int, len(r.contigs))
	var a, b int
	nIntraLinks := 0
	nInterLinks := 0
	nSkippedTooShort := 0
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
			// fmt.Println(r.contigs[ci], a, nIntraLinks, nInterLinks)
		}

		link := abs(a - b)
		if link < MinLinkDist {
			nSkippedTooShort++
			continue
		}

		// For intra-contig link it's easy, just store the distance between two ends
		// An intra-contig link
		if checkInRange(b, r.contigs[ci].start, r.contigs[ci].end) {
			r.seq.links = append(r.seq.links, link)
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
		r.interLinksFwd[ci] = append(r.interLinksFwd[ci], link)
		// Assuming flipped orientation
		link = abs(r.contigs[ci].start + r.contigs[ci].end - a - b)
		r.interLinksRev[ci] = append(r.interLinksRev[ci], link)
		nInterLinks++
	}
	log.Noticef("A total of %d intra-contig and %d inter-contig links imported (%d skipped, too short)",
		nIntraLinks, nInterLinks, nSkippedTooShort)
}

// ComputeLikelihood computes the likelihood of link sizes assuming + orientation
// and - orientation, respectively
func (r *Assesser) computeLikelihood(links []int) float64 {
	sumLogP := 0.0
	for _, link := range links {
		// if link < MinLinkDist {
		// 	link = MinLinkDist
		// }
		// bin := linkBin(link)
		sumLogP += r.model.tranformLogProb(link)
	}
	return sumLogP
}

// posterioProbability calculates the posterior probability given two likelihood
func posteriorProbability(L1, L2 float64) float64 {
	base := L1 // smaller of the two
	if L1 > L2 {
		base = L2
	}
	if L1-base > 100 {
		return 1.0
	}
	if L2-base > 100 {
		return 0.0
	}
	p1 := math.Exp(L1 - base)
	p2 := math.Exp(L2 - base)
	return p1 / (p1 + p2)
}

// computePosteriorProb computes the posterior probability of the orientations
func (r *Assesser) computePosteriorProb() {
	nFwdBetter := 0
	nRevBetter := 0
	nHighConf := 0
	r.postprob = make([]float64, len(r.contigs))
	for i := range r.contigs {
		// if i > 100 {
		// 	break
		// }
		// fmt.Println(contig)
		fwdLogP := r.computeLikelihood(r.interLinksFwd[i])
		revLogP := r.computeLikelihood(r.interLinksRev[i])
		fwdProb := posteriorProbability(fwdLogP, revLogP)
		r.postprob[i] = fwdProb
		// fmt.Println(fwdLogP, revLogP, fwdProb)

		if fwdLogP > revLogP {
			nFwdBetter++
			if fwdProb > probCutoff {
				nHighConf++
			}
		} else {
			nRevBetter++
		}
	}
	log.Noticef("Same direction better: %d; Opposite direction better: %d",
		nFwdBetter, nRevBetter)
	log.Noticef("High confidence (threshold %d%%): %d out of %d (%.1f%%)",
		int(probCutoff*100.), nHighConf, len(r.contigs), float64(nHighConf)*100./float64(len(r.contigs)))
}
