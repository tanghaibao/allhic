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
	"regexp"
	"sort"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

var b = [...]uint{0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000}
var s = [...]uint{1, 2, 4, 8, 16}

// We fit 16 bins into each power of 2
var linkRange = math.Log2(float64(MaxLinkDist) / float64(MinLinkDist))
var nBins = int(math.Ceil(linkRange * 16))

// Extracter processes the distribution step
type Extracter struct {
	Bamfile         string
	Fastafile       string
	RE              string
	MinLinks        int
	contigs         []*ContigInfo
	contigToIdx     map[string]int
	model           *LinkDensityModel
	totalIntraLinks int
	// Output file
	OutContigsfile string
	OutPairsfile   string
	OutClmfile     string
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
	label          string // allelic/cross-allelic/ok
}

// Pattern is a string pattern that is either simple or a regex
type Pattern struct {
	pattern   []byte
	rePattern *regexp.Regexp
	isRegex   bool
}

// String outputs the string representation of ContigInfo
func (r ContigInfo) String() string {
	return fmt.Sprintf("%s\t%d\t%d", r.name, r.recounts, r.length)
}

// String outputs the string representation of ContigInfo
func (r ContigPair) String() string {
	return fmt.Sprintf("%d\t%d\t%s\t%s\t%d\t%d\t%d\t%.1f\t%s",
		r.ai, r.bi, r.at, r.bt, r.RE1, r.RE2,
		r.nObservedLinks, r.nExpectedLinks, r.label)
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

// Run calls the distribution steps
func (r *Extracter) Run() {
	r.readFastaAndWriteRE()
	r.extractContigLinks()
	r.makeModel(RemoveExt(r.Bamfile) + ".distribution.txt")
	r.calcIntraContigs()
	r.calcInterContigs()
	log.Notice("Success")
}

// makeModel computes the norms and bins separately to derive an empirical link size
// distribution, then power law is inferred for extrapolating higher values
func (r *Extracter) makeModel(outfile string) {
	contigSizes := []int{}
	for _, contig := range r.contigs {
		contigSizes = append(contigSizes, contig.length)
	}
	m := NewLinkDensityModel()
	m.makeBins()
	m.makeNorms(contigSizes)
	m.countBinDensities(r.contigs)
	m.writeDistribution(outfile)
	r.model = m
}

// writeRE write a RE file and report statistics
func writeRE(outfile string, contigs []*ContigInfo) {
	f, err := os.Create(outfile)
	ErrorAbort(err)
	w := bufio.NewWriter(f)
	defer f.Close()
	totalCounts := 0
	totalBp := int64(0)
	fmt.Fprintf(w, REHeader)
	for _, contig := range contigs {
		totalCounts += contig.recounts
		totalBp += int64(contig.length)
		fmt.Fprintf(w, "%s\n", contig)
	}
	w.Flush()
	log.Noticef("RE counts in %d contigs (total: %d, avg 1 per %d bp) written to `%s`",
		len(contigs), totalCounts, totalBp/int64(totalCounts), outfile)
}

// MakePattern builds a regex-aware pattern that could be passed around and counted
// Multiple patterns will be split at comma (,) and N is converted to [ACGT]
func MakePattern(s string) Pattern {
	rePatternStr := s
	isRegex := false
	if strings.Contains(s, ",") {
		rePatternStr = ""
		for i, pattern := range strings.Split(s, ",") {
			if i != 0 {
				rePatternStr += "|"
			}
			rePatternStr += fmt.Sprintf("(%s)", pattern)
		}
		isRegex = true
	}
	if strings.Contains(s, "N") {
		rePatternStr = strings.ReplaceAll(rePatternStr, "N", "[ACGT]")
		isRegex = true
	}
	rePattern := regexp.MustCompile(rePatternStr)
	if isRegex {
		log.Noticef("Compile '%s' => '%s'", s, rePatternStr)
	}
	return Pattern{
		pattern:   []byte(s),
		rePattern: rePattern,
		isRegex:   isRegex,
	}
}

// CountPattern count how many times a pattern occurs in seq
func CountPattern(seq []byte, pattern Pattern) int {
	if pattern.isRegex {
		if all := pattern.rePattern.FindAllIndex(seq, -1); all != nil {
			return len(all)
		}
		return 0
	}
	return bytes.Count(seq, pattern.pattern)
}

// readFastaAndWriteRE writes out the number of restriction fragments, one per line
func (r *Extracter) readFastaAndWriteRE() {
	outfile := RemoveExt(r.Bamfile) + ".counts_" + strings.ReplaceAll(r.RE, ",", "_") + ".txt"
	r.OutContigsfile = outfile
	mustExist(r.Fastafile)
	reader, _ := fastx.NewDefaultReader(r.Fastafile)
	seq.ValidateSeq = false // This flag makes parsing FASTA much faster

	r.contigs = []*ContigInfo{}
	r.contigToIdx = map[string]int{}
	totalCounts := 0
	totalBp := int64(0)
	pattern := MakePattern(r.RE)

	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}

		name := string(rec.Name)
		// Strip the sequence name to get the first part up to empty space
		name = strings.Fields(name)[0]
		// Add pseudo-count of 1 to prevent division by zero
		count := CountPattern(rec.Seq.Seq, pattern) + 1
		length := rec.Seq.Length()
		totalCounts += count
		totalBp += int64(length)
		contig := &ContigInfo{
			name:     name,
			recounts: count, // To account for contigs with 0 RE sites
			length:   length,
		}

		r.contigToIdx[name] = len(r.contigs)
		r.contigs = append(r.contigs, contig)
	}
	writeRE(outfile, r.contigs)
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

// calcInterContigs calculates the MLE of distance between all contigs
func (r *Extracter) calcInterContigs() {
	clmfile := RemoveExt(r.Bamfile) + ".clm"
	lines := readClmLines(clmfile)
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
				L1: L1, L2: L2, label: "ok"}
			cp.nExpectedLinks = sumf(r.findExpectedInterContigLinks(0, L1, L2))
			cp.nObservedLinks = len(line.links)
			contigPairs[pair] = cp
		}
	}

	outfile := RemoveExt(r.Bamfile) + ".pairs.txt"
	r.OutPairsfile = outfile
	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()
	fmt.Fprintf(w, PairsFileHeader)

	allPairs := []*ContigPair{}
	for _, c := range contigPairs {
		allPairs = append(allPairs, c)
	}
	sort.Slice(allPairs, func(i, j int) bool {
		return allPairs[i].ai < allPairs[j].ai ||
			(allPairs[i].ai == allPairs[j].ai && allPairs[i].bi < allPairs[j].bi)
	})
	for _, c := range allPairs {
		if c.nObservedLinks < r.MinLinks {
			continue
		}
		fmt.Fprintln(w, c)
	}
	w.Flush()
	log.Noticef("Contig pair analyses written to `%s`", outfile)
}

// findExpectedIntraContigLinks calculates the expected number of links within a contig
func (r *Extracter) findExpectedIntraContigLinks(L int, links []int) []float64 {
	nExpectedLinks := make([]float64, nBins)
	m := r.model

	for i := 0; i < nBins; i++ {
		binStart := m.binStarts[i]
		binStop := m.binStarts[i+1]

		if binStart >= L {
			break
		}

		left := binStart
		right := min(binStop, L)
		middleX := (left + right) / 2
		middleY := L - middleX
		nObservableLinks := (right - left) * middleY

		nExpectedLinks[i] = float64(nObservableLinks) * m.linkDensity[i]
	}

	return nExpectedLinks
}

// findExpectedInterContigLinks calculates the expected number of links between two contigs
func (r *Extracter) findExpectedInterContigLinks(D, L1, L2 int) []float64 {
	if L1 > L2 {
		L1, L2 = L2, L1
	}
	nExpectedLinks := make([]float64, nBins)
	m := r.model

	for i := 0; i < nBins; i++ {
		binStart := m.binStarts[i]
		binStop := m.binStarts[i+1]

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
		nExpectedLinks[i] = float64(nObservableLinks) * m.linkDensity[i]
	}

	return nExpectedLinks
}

// extractContigLinks converts the BAM file to .clm and .ids
func (r *Extracter) extractContigLinks() {
	fh := mustOpen(r.Bamfile)
	prefix := RemoveExt(r.Bamfile)
	clmfile := prefix + ".clm"
	r.OutClmfile = clmfile

	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, err := bam.NewReader(fh, 0)
	if br == nil {
		log.Errorf("Cannot open bamfile `%s` (%s)", r.Bamfile, err)
		os.Exit(0)
	}
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
			if nLinks < r.MinLinks {
				continue
			}
			total += nLinks
			ai, bi := pair[0], pair[1]
			at, bt := r.contigs[ai].name, r.contigs[bi].name
			fmt.Fprintf(wclm, "%s%c %s%c\t%d\t%s\n",
				at, tags[i][0], bt, tags[i][1], nLinks, arrayToString(linksWithDir, " "))
		}
	}

	wclm.Flush()
	log.Noticef("Extracted %d inter-contig groups to `%s` (total = %d, maxLinks = %d, minLinks = %d)",
		len(contigPairs), clmfile, total, maxLinks, r.MinLinks)
}
