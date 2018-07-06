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
	"bytes"
	"fmt"
	"io"
	"math"
	"os"
	"sort"

	"github.com/biogo/hts/bam"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

var b = [...]uint{0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000}
var s = [...]uint{1, 2, 4, 8, 16}

// Extracter processes the distribution step
type Extracter struct {
	Bamfile         string
	Fastafile       string
	RE              string
	contigs         []*ContigInfo
	contigToIdx     map[string]int
	maxLinkDist     int
	nBins           int
	binStarts       []int
	linkDensity     []float64
	nLinks          []int
	totalIntraLinks int
	binNorms        []int
}

// ContigInfo stores results calculated from f
type ContigInfo struct {
	name           string
	recounts       int
	length         int
	links          []int // only intra-links are included in this field
	nExpectedLinks float64
	nObservedLinks int
	skip           bool
}

// ContigPair stores results calculated from findDistanceBetweenContigs
type ContigPair struct {
	ai, bi         int
	at, bt         string
	RE1, RE2       int
	L1, L2         int
	nObservedLinks int
	nExpectedLinks float64
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
	r.readFastaAndWriteRE()
	r.extractContigLinks()
	r.writeDistribution()
	r.calcIntraContigs()
	r.calcInterContigs()
}

// readFastaAndWriteRE writes out the number of restriction fragments, one per line
func (r *Extracter) readFastaAndWriteRE() {
	outfile := RemoveExt(r.Bamfile) + ".counts_" + r.RE + ".txt"

	reader, _ := fastx.NewDefaultReader(r.Fastafile)
	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()
	totalCounts := 0
	totalBp := int64(0)
	r.contigs = []*ContigInfo{}
	r.contigToIdx = map[string]int{}

	fmt.Fprintf(w, "#Contig\tRECounts\tLength\n")
	seq.ValidateSeq = false // This flag makes parsing FASTA much faster

	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}

		name := string(rec.Name)
		count := bytes.Count(rec.Seq.Seq, []byte(r.RE))
		length := rec.Seq.Length()
		totalCounts += count
		totalBp += int64(length)
		fmt.Fprintf(w, "%s\t%d\t%d\n", name, count, length)
		contig := &ContigInfo{
			name:     name,
			recounts: count, // To account for contigs with 0 RE sites
			length:   length,
		}
		r.contigToIdx[name] = len(r.contigs)
		r.contigs = append(r.contigs, contig)
	}
	w.Flush()
	log.Noticef("RE counts (total: %d, avg 1 per %d bp) written to `%s`",
		totalCounts, totalBp/int64(totalCounts), outfile)
}

// makebins makes geometric bins and count links that fall in each bin
// This heavily borrows the method form LACHESIS
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeExtracter.cc
func (r *Extracter) makebins() {
	// Step 1: make geometrically sized bins
	// We fit 16 bins into each power of 2
	linkRange := math.Log2(float64(MaxLinkDist) / float64(MinLinkDist))
	nBins := int(math.Ceil(linkRange * 16))
	r.nBins = nBins

	r.maxLinkDist = math.MinInt32
	for _, contig := range r.contigs {
		for _, link := range contig.links {
			if link > r.maxLinkDist {
				r.maxLinkDist = link
			}
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
	r.totalIntraLinks = 0
	for _, contig := range r.contigs {
		for j := 0; j < nIntraContigBins; j++ {
			z := contig.length - r.binStarts[j]
			if z < 0 {
				break
			}
			r.binNorms[j] += z
		}
		// Step 3: loop through all links and tabulate the counts
		for _, link := range contig.links {
			bin := r.LinkBin(link)
			if bin == -1 {
				continue
			}
			r.nLinks[bin]++
			r.totalIntraLinks++
		}
	}

	// Step 4: normalize to calculate link density
	r.linkDensity = make([]float64, nBins)
	for i := 0; i < nIntraContigBins; i++ {
		r.linkDensity[i] = float64(r.nLinks[i]) * BinNorm / float64(r.binNorms[i]) / float64(r.BinSize(i))
	}

	// Step 5: assume the distribution approximates 1/x for large x
	topBin := nIntraContigBins - 1
	nTopLinks := 0
	nTopLinksNeeded := r.totalIntraLinks / 100
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

// writeDistribution writes the link
func (r *Extracter) writeDistribution() {
	r.makebins()

	outfile := RemoveExt(r.Bamfile) + ".distribution.txt"
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

// String() outputs the string representation of ContigInfo
func (r ContigInfo) String() string {
	return fmt.Sprintf("%s\t%d\t%d", r.name, r.recounts, r.length)
}

// calcIntraContigs determine the local enrichment of links on this contig.
func (r *Extracter) calcIntraContigs() {
	for _, contig := range r.contigs {
		L := contig.length
		links := contig.links
		nObservedLinks := len(links)
		nExpectedLinks := r.findExpectedIntraContigLinks(L, links)
		contig.nExpectedLinks = sumf(nExpectedLinks)
		contig.nObservedLinks = nObservedLinks
	}
}

// String() outputs the string representation of ContigInfo
func (r ContigPair) String() string {
	return fmt.Sprintf("%d\t%d\t%s\t%s\t%d\t%d\t%d\t%.1f",
		r.ai, r.bi, r.at, r.bt, r.RE1, r.RE2,
		r.nObservedLinks, r.nExpectedLinks)
}

// calcInterContigs calculates the MLE of distance between all contigs
func (r *Extracter) calcInterContigs() {
	clmfile := RemoveExt(r.Bamfile) + ".clm"
	lines := ParseClmLines(clmfile)
	contigPairs := make(map[[2]int]*ContigPair)

	for i := 0; i < len(lines); i++ {
		line := &lines[i]
		at, bt := line.at, line.bt
		ai, _ := r.contigToIdx[at]
		bi, _ := r.contigToIdx[bt]

		pair := [2]int{ai, bi}
		cp, ok := contigPairs[pair]
		if !ok {
			ca, cb := r.contigs[ai], r.contigs[bi]
			L1, L2 := ca.length, cb.length
			cp = &ContigPair{ai: ai, bi: bi, at: at, bt: bt,
				RE1: ca.recounts, RE2: cb.recounts,
				L1: L1, L2: L2}
			cp.nExpectedLinks = sumf(r.findExpectedInterContigLinks(0, L1, L2))
			cp.nObservedLinks = len(line.links)
			contigPairs[pair] = cp
		}
	}

	outfile := RemoveExt(r.Bamfile) + ".pairs.txt"
	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()
	fmt.Fprintf(w, "#X\tY\tContig1\tContig2\tRE1\tRE2\tObservedLinks\tExpectedLinksIfAdjacent\n")

	allPairs := []*ContigPair{}
	for _, c := range contigPairs {
		allPairs = append(allPairs, c)
	}
	sort.Slice(allPairs, func(i, j int) bool {
		return allPairs[i].ai < allPairs[j].ai ||
			(allPairs[i].ai == allPairs[j].ai && allPairs[i].bi < allPairs[j].bi)
	})
	for _, c := range allPairs {
		fmt.Fprintln(w, c)
	}
	w.Flush()
	log.Noticef("Contig pair analyses written to `%s`", outfile)
}

// findExpectedIntraContigLinks calculates the expected number of links within a contig
func (r *Extracter) findExpectedIntraContigLinks(L int, links []int) []float64 {
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

// findExpectedInterContigLinks calculates the expected number of links between two contigs
func (r *Extracter) findExpectedInterContigLinks(D, L1, L2 int) []float64 {
	if L1 > L2 {
		L1, L2 = L2, L1
	}
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
		nExpectedLinks[i] = float64(nObservableLinks) * r.linkDensity[i] / BinNorm
	}

	return nExpectedLinks
}

// extractContigLinks converts the BAM file to .clm and .ids
func (r *Extracter) extractContigLinks() {
	fh, _ := os.Open(r.Bamfile)
	prefix := RemoveExt(r.Bamfile)
	clmfile := prefix + ".clm"

	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	fclm, _ := os.Create(clmfile)
	wclm := bufio.NewWriter(fclm)

	refs := br.Header().Refs()
	for _, ref := range refs {
		// Sanity check to see if the contig length match up between the bam and fasta
		contig := r.contigs[r.contigToIdx[ref.Name()]]
		if contig.length != ref.Len() {
			log.Errorf("Length mismatch: %s (fasta: %d bam:%d)",
				ref.Name(), contig.length, ref.Len())
		}
	}

	// Import links into pairs of contigs
	contigPairs := make(map[[2]int][][4]int)
	for {
		rec, err := br.Read()
		if err != nil {
			if err != io.EOF {
				log.Error(err)
			}
			break
		}
		// Filtering: Unmapped | Secondary | QCFail | Duplicate | Supplementary
		if rec.MapQ == 0 || rec.Flags&3844 != 0 {
			continue
		}

		// Make sure we have these contig ids
		at, bt := rec.Ref.Name(), rec.MateRef.Name()
		ai, ok := r.contigToIdx[at]
		if !ok {
			continue
		}
		bi, ok := r.contigToIdx[bt]
		if !ok {
			continue
		}

		//         read1                                               read2
		//     ---a-- X|----- dist = a2 ----|         |--- dist = b ---|X ------ b2 ------
		//     ==============================         ====================================
		//             C1 (length L1)       |----D----|         C2 (length L2)
		apos, bpos := rec.Pos, rec.MatePos
		ca, cb := r.contigs[ai], r.contigs[bi]

		// An intra-contig link
		if ai == bi {
			if link := abs(apos - bpos); link >= MinLinkDist {
				ca.links = append(ca.links, link)
			}
			continue
		}

		// An inter-contig link
		if ai > bi {
			ai, bi = bi, ai
			apos, bpos = bpos, apos
			ca, cb = cb, ca
		}

		L1 := ca.length
		L2 := cb.length
		apos2, bpos2 := L1-apos, L2-bpos
		ApBp := apos2 + bpos
		ApBm := apos2 + bpos2
		AmBp := apos + bpos
		AmBm := apos + bpos2
		pair := [2]int{ai, bi}
		contigPairs[pair] = append(contigPairs[pair], [4]int{ApBp, ApBm, AmBp, AmBm})
	}

	intraGroups := 0
	total := 0
	// Write intra-links to .dis file
	for _, contig := range r.contigs {
		if len(contig.links) == 0 {
			continue
		}
		intraGroups++
		// contig.links = unique(contig.links)
		total += len(contig.links)
	}
	log.Noticef("Extracted %d intra-contig link groups (total = %d)",
		intraGroups, total)

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
			// linksWithDir = unique(linksWithDir)
			nLinks := len(linksWithDir)
			if nLinks > maxLinks {
				maxLinks = nLinks
			}
			total += nLinks
			ai, bi := pair[0], pair[1]
			at, bt := r.contigs[ai].name, r.contigs[bi].name
			fmt.Fprintf(wclm, "%s%c %s%c\t%d\t%s\n",
				at, tags[i][0], bt, tags[i][1], nLinks, arrayToString(linksWithDir, " "))
		}
	}

	wclm.Flush()
	log.Noticef("Extracted %d inter-contig groups to `%s` (total = %d, maxLinks = %d)",
		len(contigPairs), clmfile, total, maxLinks)
}
