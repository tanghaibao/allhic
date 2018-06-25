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
	Seqsize       int
	nBins         int
	nLinks        []int
	binNorms      []int
	linkDensity   []float64
	binStarts     []int
	logP          []float64
	contigs       []BedLine
	intraLinks    []int   // Contig link sizes for all intra-contig links
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
	r.ReadBed()
	r.ExtractContigLinks()
	r.Makebins()
	r.WriteDistribution(r.Seqid + ".distribution.txt")
	r.MakeProbDist()
	r.ComputePosteriorProb()
	r.WritePostProb(r.Seqid + ".postprob.txt")
}

// WriteDistribution writes the link size distribution to file
func (r *Assesser) WriteDistribution(outfile string) {
	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()

	fmt.Fprintf(w, "#Bin\tBinStart\tBinSize\tNumLinks\tTotalSize\tLinkDensity\n")
	for i := 0; i < r.nBins; i++ {
		fmt.Fprintf(w, "%d\t%d\t%d\t%d\t%d\t%.4g\n",
			i, r.binStarts[i], r.BinSize(i), r.nLinks[i], r.binNorms[i], r.linkDensity[i])
	}

	w.Flush()
	log.Noticef("Link size distribution written to `%s`", outfile)
}

// WritePostProb writes the final posterior probability to file
func (r *Assesser) WritePostProb(outfile string) {
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

// ReadBed parses the bedfile to extract the start and stop for all the contigs
func (r *Assesser) ReadBed() {
	fh, err := os.Open(r.Bedfile)
	if err != nil {
		log.Errorf("bedfile `%s` does not exist", r.Bedfile)
		os.Exit(1)
	}
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

// ExtractContigLinks builds the probability distribution of link sizes
func (r *Assesser) ExtractContigLinks() {
	fh, _ := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	// We need the size of the SeqId to compute expected number of links
	refs := br.Header().Refs()
	for _, ref := range refs {
		if ref.Name() == r.Seqid {
			r.Seqsize = ref.Len()
		}
	}
	log.Noticef("Seq `%s` has size %d", r.Seqid, r.Seqsize)

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

		// TODO: remove this
		// if ci > 100 {
		// 	break
		// }

		link := abs(a - b)
		if link < MinLinkDist {
			nSkippedTooShort++
			continue
		}
		// For intra-contig link it's easy, just store the distance between two ends
		// An intra-contig link
		if checkInRange(b, r.contigs[ci].start, r.contigs[ci].end) {
			r.intraLinks = append(r.intraLinks, link)
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
	log.Noticef("A total of %d intra-contig links and %d inter-contig links imported (%d skipped, too short)",
		nIntraLinks, nInterLinks, nSkippedTooShort)
}

// linkBin takes a link distance and convert to a binID
func linkBin(dist int) int {
	if dist < MinLinkDist {
		return -1
	}
	distOverMin := dist / MinLinkDist
	log2i := uintLog2(uint(distOverMin))
	log2f := uintLog2Frac(float64(dist) / float64(int(MinLinkDist)<<log2i))
	return int(16*log2i + log2f)
}

// BinSize returns the size of each bin
func (r *Assesser) BinSize(i int) int {
	return r.binStarts[i+1] - r.binStarts[i]
}

// Makebins makes geometric bins and count links that fall in each bin
// This heavily borrows the method form LACHESIS
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeExtracter.cc
func (r *Assesser) Makebins() {
	// Step 1: make geometrically sized bins
	// We fit 16 bins into each power of 2
	linkRange := math.Log2(float64(MaxLinkDist) / float64(MinLinkDist))
	nBins := int(math.Ceil(linkRange * 16))

	maxLinkDist := math.MinInt32
	for _, link := range r.intraLinks {
		if link > maxLinkDist {
			maxLinkDist = link
		}
	}

	for i := 0; 16*i <= nBins; i++ {
		jpower := 1.0
		for j := 0; j < 16 && 16*i+j <= nBins; j++ {
			binStart := MinLinkDist << uint(i)
			r.binStarts = append(r.binStarts, int(float64(binStart)*jpower))
			jpower *= GeometricBinSize
		}
	}

	// Step 2: calculate assayable sequence length
	// Find the length of assayable intra-contig sequence in each bin
	intraContigLinkRange := math.Log2(float64(maxLinkDist) / float64(MinLinkDist))
	nIntraContigBins := int(math.Ceil(intraContigLinkRange * 16))

	binNorms := make([]int, nBins)
	nLinks := make([]int, nBins)
	for _, contig := range r.contigs {
		for j := 0; j < nIntraContigBins; j++ {
			z := contig.size - r.binStarts[j]
			if z < 0 {
				break
			}
			binNorms[j] += z
		}
	}

	// Step 3: loop through all links and tabulate the counts
	for _, link := range r.intraLinks {
		bin := linkBin(link)
		if bin == -1 {
			continue
		}
		nLinks[bin]++
	}

	// Step 4: normalize to calculate link density
	r.linkDensity = make([]float64, nBins)
	for i := 0; i < nIntraContigBins; i++ {
		r.linkDensity[i] = float64(nLinks[i]) / float64(binNorms[i])
	}

	// Step 5: assume the distribution approximates 1/x for large x
	topBin := nIntraContigBins - 1
	nTopLinks := 0
	nTopLinksNeeded := len(r.intraLinks) / 100
	for ; nTopLinks < nTopLinksNeeded; topBin-- {
		nTopLinks += nLinks[topBin]
	}

	avgLinkDensity := 0.0
	for i := topBin; i < nIntraContigBins; i++ {
		avgLinkDensity += r.linkDensity[i]
	}
	avgLinkDensity /= float64(nIntraContigBins - topBin)

	// Overwrite the values of last few bins, or a bin with na values
	for i := 0; i < nBins; i++ {
		if r.linkDensity[i] == 0 || i >= topBin {
			r.linkDensity[i] = avgLinkDensity
		}
	}

	// Save up for exporting to file
	r.nBins = nBins
	r.binNorms = binNorms
	r.nLinks = nLinks
}

// MakeProbDist calculates the expected number of links in each bin, which is then
// normalized to make a probability distribution
//
// TODO: check if we need to build prob dist for each contig
func (r *Assesser) MakeProbDist() {
	nBins := len(r.linkDensity)
	expectedLinks := make([]float64, nBins) // expected number of links for each bin
	sum := 0.0
	minLink := math.MaxFloat64
	for i := 0; i < nBins; i++ {
		if r.Seqsize > r.binStarts[i] {
			expectedLinks[i] = r.linkDensity[i] * float64(r.Seqsize-r.binStarts[i])
			if expectedLinks[i] < minLink {
				minLink = expectedLinks[i]
			}
		} else {
			expectedLinks[i] = minLink
		}
		sum += expectedLinks[i]
	}

	r.logP = make([]float64, nBins)
	// Normalize so that they add to 1
	for i := 0; i < nBins; i++ {
		r.logP[i] = math.Log(expectedLinks[i] / sum / float64(r.BinSize(i)))
	}
	// fmt.Println(expectedLinks)
	// fmt.Println(r.logP)
}

// ComputeLikelihood computes the likelihood of link sizes assuming + orientation
// and - orientation, respectively
func computeLikelihood(links []int, logP []float64) float64 {
	sumLogP := 0.0
	for _, link := range links {
		if link < MinLinkDist {
			link = MinLinkDist
		}
		bin := linkBin(link)
		sumLogP += logP[bin]
	}
	return sumLogP
}

// posterioProbability calculates the posterior probability given two likelihood
func posterioProbability(L1, L2 float64) float64 {
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

// ComputePosteriorProb computes the posterior probability of the orientations
func (r *Assesser) ComputePosteriorProb() {
	nFwdBetter := 0
	nRevBetter := 0
	nHighConf := 0
	r.postprob = make([]float64, len(r.contigs))
	for i := range r.contigs {
		// if i > 100 {
		// 	break
		// }
		// fmt.Println(contig)
		fwdLogP := computeLikelihood(r.interLinksFwd[i], r.logP)
		revLogP := computeLikelihood(r.interLinksRev[i], r.logP)
		fwdProb := posterioProbability(fwdLogP, revLogP)
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
