/**
 * Filename: /Users/bao/code/allhic/distribution.go
 * Path: /Users/bao/code/allhic
 * Created Date: Wednesday, March 7th 2018, 1:56:45 pm
 * Author: bao
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
	"strconv"
	"strings"

	"github.com/biogo/hts/bam"
)

var b = [...]uint{0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000}
var s = [...]uint{1, 2, 4, 8, 16}

// Distribution processes the distribution step
type Distribution struct {
	Bamfile           string
	maxLinkDist       int
	nBins             int
	binStarts         []int
	links             []int
	contigLens        []int
	linkDensity       []float64
	nLinks            []int
	binNorms          []int
	contigLinks       map[string][]int
	contigSizes       map[string]int
	contigEnrichments map[string]float64
}

// uintLog2 calculates the integer log2 of a number
func uintLog2(x2 uint) uint {
	var answer uint
	if (x2 & b[4]) != 0 {
		x2 >>= s[4]
		answer |= s[4]
	}
	if (x2 & b[3]) != 0 {
		x2 >>= s[3]
		answer |= s[3]
	}
	if (x2 & b[2]) != 0 {
		x2 >>= s[2]
		answer |= s[2]
	}
	if (x2 & b[1]) != 0 {
		x2 >>= s[1]
		answer |= s[1]
	}
	if (x2 & b[0]) != 0 {
		x2 >>= s[0]
		answer |= s[0]
	}
	return answer
}

// uintLog2Frac calculates the fractional part of log2 of a number
func uintLog2Frac(x float64) uint {
	var l uint
	x2 := x * x
	if x2 > 2 {
		x2 /= 2
		l += 8
	}
	x2 *= x2
	if x2 > 2 {
		x2 /= 2
		l += 4
	}
	x2 *= x2
	if x2 > 2 {
		x2 /= 2
		l += 2
	}
	x2 *= x2
	if x2 > 2 {
		x2 /= 2
		l++
	}
	return l
}

// LinkBin takes a link distance and convert to a binID
func (r *Distribution) LinkBin(dist int) int {
	if dist < MinLinkDist {
		return -1
	}
	distOverMin := dist / MinLinkDist
	log2i := uintLog2(uint(distOverMin))
	log2f := uintLog2Frac(float64(dist) / float64(int(MinLinkDist)<<log2i))
	return int(16*log2i + log2f)
}

// BinSize returns the size of each bin
func (r *Distribution) BinSize(i int) int {
	return r.binStarts[i+1] - r.binStarts[i]
}

// Run calls the distribution steps
func (r *Distribution) Run() {
	r.ExtractIntraContigLinks()
	r.ExtractInterContigLinks()
	r.Makebins()
	r.WriteDistribution("distribution.txt")
	r.FindEnrichmentOnContigs("enrichment.txt")
	r.FindDistanceBetweenContigs("distance.txt")
}

// Makebins makes geometric bins and count links that fall in each bin
// This heavily borrows the method form LACHESIS
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeDistribution.cc
func (r *Distribution) Makebins() {
	// Step 1: make geometrically sized bins
	// We fit 16 bins into each power of 2
	linkRange := math.Log2(float64(MaxLinkDist) / float64(MinLinkDist))
	nBins := int(math.Ceil(linkRange * 16))
	r.nBins = nBins

	r.maxLinkDist = MinInt
	for _, link := range r.links {
		if link > r.maxLinkDist {
			r.maxLinkDist = link
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
	intraContigLinkRange := math.Log2(float64(r.maxLinkDist) / float64(MinLinkDist))
	nIntraContigBins := int(math.Ceil(intraContigLinkRange * 16))

	r.binNorms = make([]int, nBins)
	r.nLinks = make([]int, nBins)
	for _, contiglen := range r.contigLens {
		for j := 0; j < nIntraContigBins; j++ {
			z := contiglen - r.binStarts[j]
			if z < 0 {
				break
			}
			r.binNorms[j] += z
		}
	}

	// Step 3: loop through all links and tabulate the counts
	for _, link := range r.links {
		bin := r.LinkBin(link)
		if bin == -1 {
			continue
		}
		r.nLinks[bin]++
	}

	// Step 4: normalize to calculate link density
	r.linkDensity = make([]float64, nBins)
	for i := 0; i < nIntraContigBins; i++ {
		r.linkDensity[i] = float64(r.nLinks[i]) * BinNorm / float64(r.binNorms[i]) / float64(r.BinSize(i))
	}

	// Step 5: assume the distribution approximates 1/x for large x
	topBin := nIntraContigBins - 1
	nTopLinks := 0
	nTopLinksNeeded := len(r.links) / 100
	for ; nTopLinks < nTopLinksNeeded; topBin-- {
		nTopLinks += r.nLinks[topBin]
	}

	avgLinkDensity := 0.0
	for i := topBin; i < nIntraContigBins; i++ {
		avgLinkDensity += r.linkDensity[i] * float64(r.BinSize(i))
	}
	avgLinkDensity /= float64(nIntraContigBins - topBin)

	// Overwrite the values of last few bins, or a bin with na values
	for i := 0; i < nBins; i++ {
		if r.linkDensity[i] == 0 || i >= topBin {
			r.linkDensity[i] = avgLinkDensity / float64(r.BinSize(i))
		}
	}
}

// WriteDistribution writes the link size distribution to file
func (r *Distribution) WriteDistribution(outfile string) {
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

// FindEnrichmentOnContigs determine the local enrichment of links on this contig.
func (r *Distribution) FindEnrichmentOnContigs(outfile string) {
	r.contigEnrichments = make(map[string]float64)
	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()

	fmt.Fprintf(w, "#Contig\tLength\tExpected\tObserved\tLDE\n")

	for contig, L := range r.contigSizes {
		links := r.contigLinks[contig]
		nObservedLinks := len(links)
		nExpectedLinks := r.FindExpectedIntraContigLinks(L, links)
		var LDE float64
		if nObservedLinks == 0 {
			LDE = 0.0
		} else {
			LDE = float64(nObservedLinks) / sumf(nExpectedLinks)
		}
		// Cap the LDE value within [0.1, 10.0]
		if LDE < 0.1 {
			LDE = 0.1
		} else if LDE > 10.0 {
			LDE = 10.0
		}
		fmt.Fprintf(w, "%s\t%d\t%.1f\t%d\t%.4f\n",
			contig, L, sumf(nExpectedLinks), nObservedLinks, LDE)

		r.contigEnrichments[contig] = LDE
	}
	w.Flush()
	log.Noticef("Link enrichments written to `%s`", outfile)
}

// FindDistanceBetweenContigs calculates the MLE of distance between all contigs
func (r *Distribution) FindDistanceBetweenContigs(outfile string) {
	clmfile := RemoveExt(r.Bamfile) + ".clm"
	lines := ParseClmLines(clmfile)
	for _, line := range lines {
		at, bt := line.at, line.bt
		L1, _ := r.contigSizes[at]
		L2, _ := r.contigSizes[bt]
		// localLDE is weighted average of LDEs of contig A and B
		lde1, _ := r.contigEnrichments[at]
		lde2, _ := r.contigEnrichments[bt]
		// Solve: lde1^L1 * lde2^L2 = x^(L1+L2)
		localLDE := math.Exp((float64(L1)*math.Log(lde1) + float64(L2)*math.Log(lde2)) / float64(L1+L2))
		fmt.Println(at, bt, L1, L2, lde1, lde2, localLDE, line.links)
		r.FindDistanceBetweenLinks(L1, L2, localLDE, line.links)
	}
}

// FindDistanceBetweenLinks calculates the most likely inter-contig distance
// Method credit to LACHESIS src code:
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeDistribution.cc
func (r *Distribution) FindDistanceBetweenLinks(L1, L2 int, LDE float64, links []int) {
	if L1 > L2 {
		L1, L2 = L2, L1
	}
	r.LogLikelihoodD(0, L1, L2, LDE, links)
}

// LogLikelihoodD calculates the log-likelihood given the distance between contigs D
// This function gets called by FindDistanceBetweenLinks
func (r *Distribution) LogLikelihoodD(D, L1, L2 int, LDE float64, links []int) {
	// Find expected number of links per bin
	nExpectedLinks := r.FindExpectedInterContigLinks(D, L1, L2, LDE)
	nObservedLinks := make([]float64, r.nBins)
	fmt.Printf("Expected: %.1f, Observed: %d\n", sumf(nExpectedLinks), len(links))

	// Find actual number of links falling into each bin
	for _, link := range links {
		bin := r.LinkBin(link + D)
		if bin == -1 {
			continue
		}
		nObservedLinks[bin]++
	}
	// fmt.Println(nObservedLinks)
}

// FindExpectedIntraContigLinks calculates the expected number of links within a contig
func (r *Distribution) FindExpectedIntraContigLinks(L int, links []int) []float64 {
	nExpectedLinks := make([]float64, r.nBins)

	for i := 0; i < r.nBins; i++ {
		binStart := r.binStarts[i]
		binStop := r.binStarts[i+1]

		if binStart >= L {
			break
		}

		left := binStart
		right := min(binStop, L)
		middleX := (left + right) / 2
		middleY := L - middleX
		nObservableLinks := (right - left) * middleY

		nExpectedLinks[i] = float64(nObservableLinks) * r.linkDensity[i] / BinNorm
	}

	return nExpectedLinks
}

// FindExpectedInterContigLinks calculates the expected number of links between two contigs
func (r *Distribution) FindExpectedInterContigLinks(D, L1, L2 int, LDE float64) []float64 {
	nExpectedLinks := make([]float64, r.nBins)

	for i := 0; i < r.nBins; i++ {
		binStart := r.binStarts[i]
		binStop := r.binStarts[i+1]

		if binStop <= D {
			continue
		}
		if binStart >= D+L1+L2 {
			break
		}

		nObservableLinks := 0

		// If the bin falls along the left slope of the trapezoid
		if binStart < D+L1 {
			left := max(binStart, D)
			right := min(binStop, D+L1)
			middleX := (left + right) / 2
			middleY := middleX - D
			nObservableLinks += (right - left) * middleY
		}

		// If the bin falls along the flat middle of the trapezoid
		if binStop >= D+L1 && binStart < D+L2 {
			left := max(binStart, D+L1)
			right := min(binStop, D+L2)
			nObservableLinks += (right - left) * L1
		}

		// If the bin falls along the right slope of the trapezoid
		if binStop >= D+L2 {
			left := max(binStart, D+L2)
			right := min(binStop, D+L1+L2)
			middleX := (left + right) / 2
			middleY := D + L1 + L2 - middleX
			nObservableLinks += (right - left) * middleY
		}
		nExpectedLinks[i] = float64(nObservableLinks) * r.linkDensity[i] / BinNorm * LDE
	}

	return nExpectedLinks
}

// ExtractIntraContigLinks populates links
func (r *Distribution) ExtractIntraContigLinks() {
	disfile := RemoveExt(r.Bamfile) + ".dis"
	log.Noticef("Parse dist file `%s`", disfile)
	r.contigLinks = make(map[string][]int)
	var contigLinks []int

	file, _ := os.Open(disfile)
	reader := bufio.NewReader(file)
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		words := strings.Split(row, "\t")
		contigLinks = []int{}
		for _, link := range strings.Split(words[1], ",") {
			ll, _ := strconv.Atoi(link)
			if ll >= MinLinkDist {
				r.links = append(r.links, ll)
				contigLinks = append(contigLinks, ll)
			}
		}
		r.contigLinks[words[0]] = contigLinks
	}
	log.Noticef("Imported %d intra-contig links", len(r.links))
}

// ExtractInterContigLinks converts the BAM file to .clm and .ids
func (r *Distribution) ExtractInterContigLinks() {
	fh, _ := os.Open(r.Bamfile)
	prefix := RemoveExt(r.Bamfile)
	disfile := prefix + ".dis"
	clmfile := prefix + ".clm"
	idsfile := prefix + ".ids"

	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	fdis, _ := os.Create(disfile)
	wdis := bufio.NewWriter(fdis)
	fclm, _ := os.Create(clmfile)
	wclm := bufio.NewWriter(fclm)
	fids, _ := os.Create(idsfile)
	wids := bufio.NewWriter(fids)

	r.contigSizes = make(map[string]int)
	refs := br.Header().Refs()
	for _, ref := range refs {
		r.contigLens = append(r.contigLens, ref.Len())
		r.contigSizes[ref.Name()] = ref.Len()
		fmt.Fprintf(wids, "%s\t%d\n", ref.Name(), ref.Len())
	}
	wids.Flush()
	log.Noticef("Extracted %d contigs to `%s`", len(r.contigLens), idsfile)

	// Import links into pairs of contigs
	var a, a2, b, b2 int
	contigPairs := make(map[[2]string][][4]int)
	contigLinks := make(map[string][]int)
	for {
		rec, err := br.Read()
		if err != nil {
			if err != io.EOF {
				log.Error(err)
			}
			break
		}
		at, bt := rec.Ref.Name(), rec.MateRef.Name()

		//         read1                                               read2
		//     ---a-- X|----- dist = a2 ----|         |--- dist = b ---|X ------ b2 ------
		//     ==============================         ====================================
		//             C1 (length L1)       |----D----|         C2 (length L2)
		rlen := rec.Len()
		a, b = rec.Pos, rec.MatePos

		// An intra-contig link
		if at == bt {
			if link := abs(a - b); link >= MinLinkDist {
				contigLinks[at] = append(contigLinks[at], link)
			}
			continue
		}

		// An inter-contig link
		if at > bt {
			at, bt = bt, at
			a, b = b, a
		}

		L1, _ := r.contigSizes[at]
		L2, _ := r.contigSizes[bt]
		a2, b2 = L1-rlen-a, L2-rlen-b
		ApBp := a2 + b
		ApBm := a2 + b2
		AmBp := a + b
		AmBm := a + b2
		pair := [2]string{at, bt}
		contigPairs[pair] = append(contigPairs[pair], [4]int{ApBp, ApBm, AmBp, AmBm})
	}

	// Write intra-links to .dis file
	for contig, links := range contigLinks {
		links = unique(links)
		fmt.Fprintf(wdis, "%s\t%s\n", contig, arrayToString(links, ","))
	}
	wdis.Flush()
	log.Noticef("Extracted %d intra-contig links to `%s`", len(contigLinks), disfile)

	// Write inter-links to .clm file
	tags := []string{"++", "+-", "-+", "--"}
	for pair, links := range contigPairs {
		for i := 0; i < 4; i++ {
			linksWithDir := make([]int, len(links))
			for j, link := range links {
				linksWithDir[j] = link[i]
			}
			linksWithDir = unique(linksWithDir)
			fmt.Fprintf(wclm, "%s%c %s%c\t%d\t%s\n",
				pair[0], tags[i][0], pair[1], tags[i][1], len(linksWithDir), arrayToString(linksWithDir, " "))
		}
	}

	wclm.Flush()
	log.Noticef("Extracted %d inter-contig links to `%s`", len(contigPairs), clmfile)
}
