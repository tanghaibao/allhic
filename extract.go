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

// Extracter processes the distribution step
type Extracter struct {
	Bamfile           string
	maxLinkDist       int
	nBins             int
	logFactorials     []float64
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
func (r *Extracter) LinkBin(dist int) int {
	if dist < MinLinkDist {
		return -1
	}
	distOverMin := dist / MinLinkDist
	log2i := uintLog2(uint(distOverMin))
	log2f := uintLog2Frac(float64(dist) / float64(int(MinLinkDist)<<log2i))
	return int(16*log2i + log2f)
}

// BinSize returns the size of each bin
func (r *Extracter) BinSize(i int) int {
	return r.binStarts[i+1] - r.binStarts[i]
}

// Run calls the distribution steps
func (r *Extracter) Run() {
	r.ExtractInterContigLinks()
	r.ExtractIntraContigLinks()
	r.Makebins()
	r.WriteExtracter("distribution.txt")
	r.FindEnrichmentOnContigs("enrichment.txt")
	r.PreComputeLogFactorials()
	r.FindDistanceBetweenContigs("distance.txt")
}

// Makebins makes geometric bins and count links that fall in each bin
// This heavily borrows the method form LACHESIS
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeExtracter.cc
func (r *Extracter) Makebins() {
	// Step 1: make geometrically sized bins
	// We fit 16 bins into each power of 2
	linkRange := math.Log2(float64(MaxLinkDist) / float64(MinLinkDist))
	nBins := int(math.Ceil(linkRange * 16))
	r.nBins = nBins

	r.maxLinkDist = math.MinInt32
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

// WriteExtracter writes the link size distribution to file
func (r *Extracter) WriteExtracter(outfile string) {
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

// ContigInfo stores results calculated from FindEnrichmentOnContigs
type ContigInfo struct {
	name           string
	length         int
	nExpectedLinks float64
	nObservedLinks int
	lde            float64
}

// String() outputs the string representation of ContigInfo
func (r ContigInfo) String() string {
	return fmt.Sprintf("%s\t%d\t%.1f\t%d\t%.4f",
		r.name, r.length, r.nExpectedLinks, r.nObservedLinks, r.lde)
}

// FindEnrichmentOnContigs determine the local enrichment of links on this contig.
func (r *Extracter) FindEnrichmentOnContigs(outfile string) {
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
		// Cap the LDE value within [0.2, 5.0]
		if LDE < 0.2 {
			LDE = 0.2
		} else if LDE > 5.0 {
			LDE = 5.0
		}
		ci := ContigInfo{name: contig, length: L, nExpectedLinks: sumf(nExpectedLinks),
			nObservedLinks: nObservedLinks, lde: LDE}
		fmt.Fprintln(w, ci)

		r.contigEnrichments[contig] = LDE
	}
	w.Flush()
	log.Noticef("Link enrichments written to `%s`", outfile)
}

// ContigPair stores results calculated from FindDistanceBetweenContigs
type ContigPair struct {
	at             string
	bt             string
	L1             int
	L2             int
	lde1           float64
	lde2           float64
	localLDE       float64
	nObservedLinks int
	nExpectedLinks float64
	mleDistance    int
	logLikelihood  float64
	score          float64
}

// String() outputs the string representation of ContigInfo
func (r ContigPair) String() string {
	return fmt.Sprintf("%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\t%.1f\t%d\t%.1f",
		r.at, r.bt, r.L1, r.L2, r.lde1, r.lde2, r.localLDE,
		r.nObservedLinks, r.nExpectedLinks, r.mleDistance, r.score)
}

// FindDistanceBetweenContigs calculates the MLE of distance between all contigs
func (r *Extracter) FindDistanceBetweenContigs(outfile string) {
	clmfile := RemoveExt(r.Bamfile) + ".clm"
	lines := ParseClmLines(clmfile)
	contigPairs := make(map[[2]string]*ContigPair)
	var L1, L2 int
	var lde1, lde2, localLDE float64

	longestContigSize := 0
	for _, v := range r.contigSizes {
		if v > longestContigSize {
			longestContigSize = v
		}
	}
	longestContigSizeSquared := float64(longestContigSize) * float64(longestContigSize)

	for i := 0; i < len(lines); i++ {
		line := &lines[i]
		at, bt := line.at, line.bt
		pair := [2]string{at, bt}
		cp, cpok := contigPairs[pair]
		if !cpok {
			L1, _ = r.contigSizes[at]
			L2, _ = r.contigSizes[bt]
			// localLDE is weighted average of LDEs of contig A and B
			lde1, _ = r.contigEnrichments[at]
			lde2, _ = r.contigEnrichments[bt]
			// Solve: lde1^L1 * lde2^L2 = x^(L1+L2)
			localLDE = math.Exp((float64(L1)*math.Log(lde1) + float64(L2)*math.Log(lde2)) / float64(L1+L2))
			cp = &ContigPair{at: at, bt: bt, L1: L1, L2: L2, lde1: lde1, lde2: lde2, localLDE: localLDE}
			contigPairs[pair] = cp
			r.FindDistanceBetweenLinks(cp, line)
		}
	}

	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()
	fmt.Fprintf(w, "#Contig1\tContig2\tLength1\tLength2\tLDE1\tLDE2\tLDE"+
		"\tObservedLinks\tExpectedLinksIfAdjacent\tMLEdistance\tNormalizedScore\n")
	for _, c := range contigPairs {
		c.score = float64(c.nObservedLinks) * longestContigSizeSquared / (float64(c.L1) * float64(c.L2))
		fmt.Fprintln(w, c)
	}
	w.Flush()
	log.Noticef("Contig pair analyses written to `%s`", outfile)
}

// FindDistanceBetweenLinks calculates the most likely inter-contig distance
// Method credit to LACHESIS src code:
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeExtracter.cc
func (r *Extracter) FindDistanceBetweenLinks(cp *ContigPair, line *CLMLine) {
	L1, L2 := cp.L1, cp.L2
	LDE := cp.localLDE
	links := line.links

	if L1 > L2 {
		L1, L2 = L2, L1
	}

	bestLogLikelihood := math.Inf(-1)
	bestD := MaxLinkDist
	for i := 0; i < r.nBins; i++ {
		D := r.binStarts[i]
		logLikelihood := r.LogLikelihoodD(D, L1, L2, LDE, links)
		if logLikelihood > bestLogLikelihood {
			bestLogLikelihood = logLikelihood
			bestD = D
		}
	}

	if cp.logLikelihood == 0 || cp.logLikelihood < bestLogLikelihood {
		cp.mleDistance, cp.logLikelihood = bestD, bestLogLikelihood
		cp.nExpectedLinks = sumf(r.FindExpectedInterContigLinks(0, L1, L2, LDE))
		cp.nObservedLinks = len(links)
	}
}

// LogLikelihoodD calculates the log-likelihood given the distance between contigs D
// This function gets called by FindDistanceBetweenLinks
func (r *Extracter) LogLikelihoodD(D, L1, L2 int, LDE float64, links []int) float64 {
	// Find expected number of links per bin
	nExpectedLinks := r.FindExpectedInterContigLinks(D, L1, L2, LDE)
	// nObservedLinks := make([]int, r.nBins)
	// // Find actual number of links falling into each bin
	// for _, link := range links {
	// 	bin := r.LinkBin(link + D)
	// 	if bin < 0 || bin >= r.nBins {
	// 		continue
	// 	}
	// 	nObservedLinks[bin]++
	// }

	// Calculate the log likelihood of the observed link distribution
	// The log likelihood of the data is sum of log likelihood of all bins where
	// logLL = -m + k ln m - ln(k!)
	// logLikelihood := 0.0
	// for i := 0; i < r.nBins; i++ {
	// 	m := nExpectedLinks[i]
	// 	k := nObservedLinks[i]

	// 	if m == 0 { // empty bin
	// 		continue
	// 	}
	// 	logLikelihood += -m + float64(k)*math.Log(float64(m)) - r.logFactorials[k]
	// }
	m := sumf(nExpectedLinks)
	k := len(links)
	logLikelihood := -m + float64(k)*math.Log(float64(m)) - r.logFactorials[k]

	return logLikelihood
}

// FindExpectedIntraContigLinks calculates the expected number of links within a contig
func (r *Extracter) FindExpectedIntraContigLinks(L int, links []int) []float64 {
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
func (r *Extracter) FindExpectedInterContigLinks(D, L1, L2 int, LDE float64) []float64 {
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
func (r *Extracter) ExtractIntraContigLinks() {
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
func (r *Extracter) ExtractInterContigLinks() {
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
	var a, a2, b, b2, total int
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
		total += len(links)
		fmt.Fprintf(wdis, "%s\t%s\n", contig, arrayToString(links, ","))
	}
	wdis.Flush()
	log.Noticef("Extracted %d intra-contig link groups to `%s` (total = %d)",
		len(contigLinks), disfile, total)

	// Write inter-links to .clm file
	total = 0
	maxLinks := 0
	tags := []string{"++", "+-", "-+", "--"}
	for pair, links := range contigPairs {
		for i := 0; i < 4; i++ {
			linksWithDir := make([]int, len(links))
			for j, link := range links {
				linksWithDir[j] = link[i]
			}
			linksWithDir = unique(linksWithDir)
			nLinks := len(linksWithDir)
			if nLinks < MinInterLinks {
				continue
			}
			if nLinks > maxLinks {
				maxLinks = nLinks
			}
			total += nLinks
			fmt.Fprintf(wclm, "%s%c %s%c\t%d\t%s\n",
				pair[0], tags[i][0], pair[1], tags[i][1], nLinks, arrayToString(linksWithDir, " "))
		}
	}

	wclm.Flush()
	log.Noticef("Extracted %d inter-contig groups to `%s` (total = %d, maxLinks = %d)",
		len(contigPairs), clmfile, total, maxLinks)

	r.logFactorials = make([]float64, maxLinks+1)
}

// PreComputeLogFactorials precomputes log factorials for calculation of log likelihood
func (r *Extracter) PreComputeLogFactorials() {
	for i := 2; i < len(r.logFactorials); i++ {
		r.logFactorials[i] = r.logFactorials[i-1] + math.Log(float64(i))
	}
}
