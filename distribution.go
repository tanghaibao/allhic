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
	Distfile          string
	Clmfile           string
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
	r.ExtractLinks()
	r.ExtractContigLens()
	r.Makebins()
	r.WriteDistribution("distribution.txt")
	r.FindEnrichmentOnContigs("enrichment.txt")
}

// ExtractLinks populates links
func (r *Distribution) ExtractLinks() {
	file, _ := os.Open(r.Distfile)
	log.Noticef("Parse dist file `%s`", r.Distfile)
	r.contigLinks = make(map[string][]int)
	var contigLinks []int

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

// ExtractContigLens populates contigLens
func (r *Distribution) ExtractContigLens() {
	fh, _ := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	r.contigSizes = make(map[string]int)
	refs := br.Header().Refs()
	for _, ref := range refs {
		r.contigLens = append(r.contigLens, ref.Len())
		r.contigSizes[ref.Name()] = ref.Len()
	}
	log.Noticef("Imported %d contigs", len(r.contigLens))
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
	// var contigLens = [...]int{43270923, 35937250, 36413819, 35502694, 29958434, 31248787,
	// 	29697621, 28443022, 23012720, 23207287, 29021106, 27531856, 633585, 592136}

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

	// Step 6: Serialize the distribution
	sum := 0.0
	for i, ld := range r.linkDensity {
		sum += ld * float64(r.BinSize(i))
	}

	normLinkDensity := make([]float64, nBins)
	copy(normLinkDensity, r.linkDensity)
	for i, ld := range r.linkDensity {
		normLinkDensity[i] = ld * float64(r.BinSize(i)) / sum
	}
	fmt.Println(normLinkDensity)
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
		LDE := float64(nObservedLinks) / sumf(nExpectedLinks)
		fmt.Fprintf(w, "%s\t%d\t%.1f\t%d\t%.1f\n",
			contig, L, sumf(nExpectedLinks), nObservedLinks, LDE)

		r.contigEnrichments[contig] = LDE
	}
	w.Flush()
	log.Noticef("Link enrichments written to `%s`", outfile)
}

// FindDistanceBetweenLinks calculates the most likely inter-contig distance
// Method credit to LACHESIS src code:
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeDistribution.cc
func (r *Distribution) FindDistanceBetweenLinks(sizeA, sizeB int, LDE float64, links []int) {

}

// LogLikelihoodD calculates the log-likelihood given the distance between contigs D
// This function gets called by FindDistanceBetweenLinks
func (r *Distribution) LogLikelihoodD(D, L1, L2 int, LDE float64, links []int) {
	// Step 1. Find expected number of links per bin
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
func (r *Distribution) FindExpectedInterContigLinks(D, L1, L2 int, LDE float64, links []int) []float64 {
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
	}

	return nExpectedLinks
}
