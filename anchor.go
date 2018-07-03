/*
 * Filename: /Users/bao/code/allhic/anchor.go
 * Path: /Users/bao/code/allhic
 * Created Date: Monday, June 4th 2018, 9:26:26 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/kshedden/gonpy"
)

// Anchorer runs the merging algorithm
type Anchorer struct {
	Bamfile      string
	contigs      []*Contig
	nameToContig map[string]*Contig
	path         *Path
}

// AnchorerJSON keeps a succinct subset of all fields in Anchorer
type AnchorerJSON struct {
	Starts        map[string]int64 `json:"starts"`
	Sizes         map[string]int64 `json:"sizes"`
	TotalBins     int              `json:"total_bins"`
	DistBinStarts []int64          `json:"distbinstarts"`
	DistBinSizes  []int64          `json:"distbinsizes"`
	Resolution    int64            `json:"resolution"`
}

// Contig stores the name and length of each contig
type Contig struct {
	name        string
	length      int64
	links       []*Link
	path        *Path
	start       int64
	orientation int8 // 1 => +, -1 => -
	segments    []Range
}

// Link contains a specific inter-contig link
type Link struct {
	a, b       *Contig // Contig ids
	apos, bpos int64   // Positions
}

// Path is a collection of ordered contigs
type Path struct {
	contigs      []*Contig // List of contigs
	LNode, RNode *Node     // Two nodes at each end
	length       int64     // Cumulative length of all contigs
}

// Range tracks contig:start-end
type Range struct {
	start int64
	end   int64
	node  *Node
}

// Piler has the data structures to support overlap counts
// Here we use a data structure described in:
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530906/
// We store the starts and ends of links in sorted arrays
// The `icount` algorithm then search an interval (or in this case) point
// query into these sorted interval ends
type Piler struct {
	BS, BE []int64
}

// SparseMatrix stores a big square matrix that is sparse
type SparseMatrix []map[int]int

// PathSet stores the set of paths
type PathSet map[*Path]bool

// Run kicks off the merging algorithm
func (r *Anchorer) Run() {
	// Prepare the paths to run
	nIterations := 1
	r.ExtractInterContigLinks()
	flanksize := int64(1000000)
	paths := r.makeTrivialPaths(r.contigs, flanksize)
	for i := 0; i < nIterations; i++ {
		r.iterativeGraphMerge(paths, flanksize)
		// Attempt to split bad joins
		// r.makeContigStarts()
		// piler := r.inspectGaps(LinkDist)
		// countCutoff := r.findCountCutoff(piler)
		// paths = r.splitPath(piler, countCutoff, flanksize)
	}

	// Serialize to disk for plotting
	r.makeContigStarts()
	r.serialize(250000, "genome.json", "data.npy")

	// printPaths(paths)
	log.Notice("Success")
}

func (r *Anchorer) iterativeGraphMerge(paths PathSet, flanksize int64) {
	var G Graph
	i := 0
	prevPaths := len(paths)
	graphRemake := true
	for prevPaths > 20 {
		// flanksize = getL50(paths) / 4
		if graphRemake {
			i++
			log.Noticef("Starting iteration %d with %d paths (L50=%d)",
				i, len(paths), getL50(paths))
			G = r.makeGraph(paths)
		}
		CG := r.makeConfidenceGraph(G)
		paths = r.generatePathAndCycle(CG, flanksize)
		// Check if no merges were made in this round
		if len(paths) == prevPaths {
			paths = r.removeSmallestPath(paths, G)
			graphRemake = false
		} else {
			graphRemake = true
		}
		prevPaths = len(paths)
	}

	printPaths(paths)

	// Path found
	r.path = nil
	for path := range paths {
		if r.path == nil || path.length > r.path.length {
			r.path = path
		}
	}
}

// removeSmallestPath forces removal of the smallest path so that we can continue
// with the merging. This is also a small operation so we'll have to just modify
// the graph only slightly
func (r *Anchorer) removeSmallestPath(paths PathSet, G Graph) PathSet {
	var smallestPath *Path
	for path := range paths {
		if smallestPath == nil || path.length < smallestPath.length {
			smallestPath = path
		}
	}
	// Inactivate the nodes
	log.Noticef("Inactivate path %s (length=%d)", smallestPath, smallestPath.length)

	// Un-assign the contigs
	for _, contig := range smallestPath.contigs {
		contig.path = nil
	}

	for _, node := range []*Node{smallestPath.LNode, smallestPath.RNode} {
		if nb, ok := G[node]; ok {
			for b := range nb {
				delete(G[b], node)
			}
			delete(G, node)
		}
	}
	delete(paths, smallestPath)
	return paths
}

// printPaths shows the current details of the clustering
func printPaths(paths PathSet) {
	for path := range paths {
		fmt.Println(path)
	}
}

// makeTrivialPaths starts the initial construction of Path object, with one
// contig per Path (trivial Path)
func (r *Anchorer) makeTrivialPaths(contigs []*Contig, flanksize int64) PathSet {
	// Initially make every contig a single Path object
	paths := PathSet{}
	for _, contig := range contigs {
		contig.orientation = 1
		makePath([]*Contig{contig}, paths, flanksize)
	}

	return paths
}

// ExtractInterContigLinks extracts links from the Bamfile
func (r *Anchorer) ExtractInterContigLinks() {
	fh, _ := os.Open(r.Bamfile)
	prefix := RemoveExt(r.Bamfile)
	disfile := prefix + ".dis"
	idsfile := prefix + ".ids"

	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	fdis, _ := os.Create(disfile)
	wdis := bufio.NewWriter(fdis)
	fids, _ := os.Create(idsfile)
	wids := bufio.NewWriter(fids)

	r.nameToContig = make(map[string]*Contig)
	refs := br.Header().Refs()
	for _, ref := range refs {
		contig := Contig{
			name:   ref.Name(),
			length: int64(ref.Len()),
		}
		r.contigs = append(r.contigs, &contig)
		r.nameToContig[contig.name] = &contig
		fmt.Fprintf(wids, "%s\t%d\n", ref.Name(), ref.Len())
	}
	wids.Flush()
	log.Noticef("Extracted %d contigs to `%s`", len(r.contigs), idsfile)

	// Import links into pairs of contigs
	intraTotal, interTotal := 0, 0
	intraLinks := make(map[string][]int)
	for {
		rec, err := br.Read()
		if err != nil {
			if err != io.EOF {
				log.Error(err)
			}
			break
		}

		at, bt := rec.Ref.Name(), rec.MateRef.Name()
		a, b := r.nameToContig[at], r.nameToContig[bt]
		apos, bpos := rec.Pos, rec.MatePos

		// An intra-contig link
		if a == b {
			if link := abs(apos - bpos); link >= MinLinkDist {
				intraLinks[at] = append(intraLinks[at], link)
			}
			intraTotal++
			continue
		}

		// An inter-contig link
		a.links = append(a.links, &Link{
			a: a, b: b, apos: int64(apos), bpos: int64(bpos),
		})
		interTotal++
	}

	for _, contig := range r.contigs {
		sort.Slice(contig.links, func(i, j int) bool {
			return contig.links[i].apos < contig.links[j].apos
		})
	}

	// Write intra-links to .dis file
	for contig, links := range intraLinks {
		links = unique(links)
		fmt.Fprintf(wdis, "%s\t%s\n", contig, arrayToString(links, ","))
	}
	wdis.Flush()
	log.Noticef("Extracted %d intra-contig and %d inter-contig links",
		intraTotal, interTotal)
}

// reverse reverses the orientations of all components
func (r *Path) reverse() {
	c := r.contigs
	for i, j := 0, len(c)-1; i < j; i, j = i+1, j-1 {
		c[i], c[j] = c[j], c[i]
	}
	for _, contig := range c {
		contig.orientation = -contig.orientation
	}
}

// String prints the Path nicely
func (r *Path) String() string {
	tagContigs := make([]string, len(r.contigs))
	for i, contig := range r.contigs {
		tag := ""
		if contig.orientation < 0 {
			tag = "-"
		}
		tagContigs[i] = tag + contig.name
	}
	return strings.Join(tagContigs, " ")
}

// contigToNode takes as input contig and position, returns the nodeID
func contigToNode(contig *Contig, pos int64) *Node {
	for _, rr := range contig.segments { // multiple 'segments'
		if pos >= rr.start && pos < rr.end {
			return rr.node
		}
	}
	// log.Errorf("%s:%d not found", contig.name, pos)
	return nil
}

// linkToNodes takes as input link, returns two nodeIDs
func (r *Anchorer) linkToNodes(link *Link) (*Node, *Node) {
	a := contigToNode(link.a, link.apos)
	b := contigToNode(link.b, link.bpos)
	return a, b
}

// insertEdge adds just one link to the graph
func (r *Anchorer) insertEdge(G Graph, a, b *Node) {
	if _, aok := G[a]; aok {
		G[a][b]++
	} else {
		G[a] = map[*Node]int64{b: 1}
	}
}

// bisect cuts the Path into two parts
// The initial implementation cuts the path into two equal halves. However, this is
// not desired for longer path since the 'internal' links should be penalized somehow.
// Therefore, we add a flanksize parameter so that we only create end nodes that are
// min(pathlength / 2, flanksize) so that we handle long paths properly
func (r *Path) bisect(flanksize int64) {
	var contig *Contig
	var contigpos int64

	r.length = 0
	for _, contig := range r.contigs {
		r.length += contig.length
	}

	flanksize = minInt64(r.length/2, flanksize)

	LNode := &Node{
		path:   r,
		length: flanksize,
	}
	RNode := &Node{
		path:   r,
		length: flanksize,
	}
	LNode.sister = RNode
	RNode.sister = LNode
	r.LNode = LNode
	r.RNode = RNode

	// Build the left node
	contigStart := int64(0)
	i := 0
	for i = 0; i < len(r.contigs); i++ {
		contig = r.contigs[i]
		if contigStart+contig.length > flanksize {
			break
		}
		contig.segments = []Range{
			Range{0, contig.length, LNode},
		}
		contigStart += contig.length
	}

	// Left flank cuts through this contig
	contigpos = flanksize - contigStart
	if contig.orientation > 0 {
		contig.segments = append(contig.segments, Range{
			0, contigpos, LNode,
		})
	} else {
		contigpos = contig.length - contigpos
		contig.segments = append(contig.segments, Range{
			contigpos, contig.length, LNode,
		})
	}

	// Fast forward to right flank
	for ; i < len(r.contigs); i++ {
		contig = r.contigs[i]
		if contigStart+contig.length > r.length-flanksize {
			break
		}
		contig.segments = nil
		contigStart += contig.length
	}

	// Right flank cuts through this contig
	contigpos = r.length - flanksize - contigStart
	if contig.orientation > 0 {
		contig.segments = append(contig.segments, Range{
			contigpos, contig.length, RNode,
		})
	} else {
		contigpos = contig.length - contigpos
		contig.segments = append(contig.segments, Range{
			0, contigpos, RNode,
		})
	}

	// Build the right node
	for ; i < len(r.contigs); i++ {
		contig = r.contigs[i]
		contig.segments = []Range{
			Range{0, contig.length, RNode},
		}
	}
}

// makeContigStarts returns starts of contigs within a path
func (r *Anchorer) makeContigStarts() {
	pos := int64(0)
	for _, contig := range r.path.contigs {
		contig.start = pos
		pos += contig.length
	}
}

// findPosition returns the i-th bin along the path
func findPosition(contig *Contig, pos int64) int64 {
	offset := pos
	if contig.orientation < 0 {
		offset = contig.length - pos
	}
	return contig.start + offset
}

// inspectGaps check each gap for the number of links <= 1Mb going across
func (r *Anchorer) inspectGaps(cutoff int64) Piler {
	// We need to quickly map all links to their [start, end] on the path
	// then increment all the link counts for each of the intervening gaps
	p := Piler{}
	for _, contig := range r.path.contigs {
		for _, link := range contig.links {
			S := findPosition(link.a, link.apos)
			E := findPosition(link.b, link.bpos)
			if S > E {
				S, E = E, S
			}
			if E-S < cutoff {
				p.BS = append(p.BS, S)
				p.BE = append(p.BE, E)
			}
		}
	}
	sortInt64s(p.BS)
	sortInt64s(p.BE)
	log.Noticef("Sorted both ends of %d links (only <%d bp included)",
		len(p.BS), LinkDist)

	return p
}

// findCountCutoff samples evenly in the genome and builds a distribution
// of the number of links going across at each position, then take the 10-th
// percentile as the cutoff
func (r *Anchorer) findCountCutoff(piler Piler) int {
	// Sample evenly in the genome to compute a cutoff value
	stepSize := r.path.length / 1000
	counts := []int{}
	for i := int64(0); i < r.path.length; i += stepSize {
		counts = append(counts, piler.intervalCounts(i))
	}
	sort.Ints(counts)
	// fmt.Println(counts)
	countCutoff := counts[len(counts)/10]
	log.Noticef("Weak joins defined as <%d links (only <%d bp included)",
		countCutoff, LinkDist)

	return countCutoff
}

// intervalCounts returns the number of intervals over the query position
func (r *Piler) intervalCounts(pos int64) int {
	return searchInt64(r.BS, pos) - searchInt64(r.BE, pos)
	// first := searchInt64(r.BS, pos)
	// last := searchInt64(r.BE, pos)
	// if first-last < 5 {
	// 	for i := max(first-5, 0); i < min(first+5, len(r.BS)); i++ {
	// 		fmt.Println(i, r.BS[i], r.BE[i])
	// 	}
	// }
	// return first - last
}

// identifyGap prints out all the gaps that lie within the bin
func (r *Anchorer) identifyGap(breakPoints []int, res int64) {
	contigStart := int64(0)
	j := 0
	for i, contig := range r.path.contigs {
		if contigStart >= int64(breakPoints[j])*res {
			if contigStart >= int64(breakPoints[j]+1)*res {
				for j < len(breakPoints) && contigStart >= int64(breakPoints[j]+1)*res {
					j++
				}
				// Exhausted, terminate
				if j == len(breakPoints) {
					break
				}
			} else {
				// We have a candidate
				fmt.Println(breakPoints[j], res, i, contigStart, contig.name)
				// fmt.Println(path.contigs[max(0, i-5):min(len(path.contigs)-1, i+5)])
			}
		}
		contigStart += contig.length
	}
}

// splitPath takes a path and look at joins that are weak
func (r *Anchorer) splitPath(piler Piler, countCutoff int, flanksize int64) PathSet {
	contigs := []*Contig{}
	paths := PathSet{}
	strength := 0
	// Now go through all the contig joins
	for _, contig := range r.path.contigs {
		// Contigs at the end of the chromosome are at a disadvantage
		if contig.start >= LinkDist && contig.start < r.path.length-LinkDist {
			// Query this join
			strength = piler.intervalCounts(contig.start)
			if strength < countCutoff { // needs to break a join here
				fmt.Println("-------------------")
				path := makePath(contigs, paths, flanksize)
				fmt.Println(path, len(path.contigs), path.length)
				contigs = []*Contig{}
			}
		}
		contigs = append(contigs, contig)
		// fmt.Println(contig.name, contig.start, contig.orientation, strength)
	}
	// Last piece
	makePath(contigs, paths, flanksize)
	// fmt.Println(path, len(path.contigs), path.length)
	log.Noticef("Split into %d paths", len(paths))

	return paths
}

// makePath creates a Path from contigs and set everything properly
func makePath(contigs []*Contig, paths PathSet, flanksize int64) *Path {
	path := &Path{
		contigs: contigs,
	}
	for _, contig := range contigs {
		contig.path = path
	}
	path.bisect(flanksize)
	paths[path] = true

	return path
}

// findBin returns the i-th bin along the path
func findBin(contig *Contig, pos, resolution int64) int {
	position := findPosition(contig, pos)
	return int(position / resolution)
}

// serialize outputs the current path to disk
// This contains the data for jcvi.assembly.hic.heatmap()
func (r *Anchorer) serialize(res int64, jsonfile, npyfile string) {
	A := &AnchorerJSON{Resolution: res}
	m := int(math.Ceil(float64(r.path.length) / float64(res)))
	A.Starts = make(map[string]int64)
	A.Sizes = make(map[string]int64)
	A.TotalBins = m
	// Initialize the count matrix
	C := make([]int32, m*m)
	for _, contig := range r.path.contigs {
		A.Starts[contig.name] = contig.start / res
		A.Sizes[contig.name] = contig.length / res
		for _, link := range contig.links {
			// The link need to be within the path!
			if link.a.path != link.b.path {
				continue
			}
			a := findBin(link.a, link.apos, res)
			b := findBin(link.b, link.bpos, res)
			C[a*m+b]++
			C[b*m+a]++
		}
	}

	// Serialize the contig size stats to JSON file
	s, _ := json.MarshalIndent(A, "", "\t")
	f, _ := os.Create(jsonfile)
	defer f.Close()
	jw := bufio.NewWriter(f)
	jw.WriteString(string(s))
	jw.Flush()
	log.Noticef("Contig stats (N=%d Length=%d) written to `%s`",
		m, r.path.length, jsonfile)

	// Serialize the pixelated matrix to NPY file
	w, _ := gonpy.NewFileWriter(npyfile)
	w.Shape = []int{m, m}
	w.WriteInt32(C)
	log.Noticef("Matrix (resolution=%d) written to `%s`", res, npyfile)
}

// getL50 computes the L50 of all component contigs within a path
func getL50(paths PathSet) int64 {
	pathLengths := []int64{}
	for path := range paths {
		pathLengths = append(pathLengths, path.length)
	}

	return L50(pathLengths)
}
