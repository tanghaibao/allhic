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
	Tourfile     string
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
func (r *Anchorer) Run() error {
	// Prepare the paths to run
	nIterations := 1
	err := r.ExtractInterContigLinks()
	if err != nil {
		return err
	}
	flanksize := int64(LIMIT)
	paths, err := r.makeTrivialPaths(r.contigs, flanksize)
	if err != nil {
		return err
	}
	for i := 0; i < nIterations; i++ {
		err := r.iterativeGraphMerge(paths, flanksize)
		if err != nil {
			return err
		}
	}

	// Serialize to disk for plotting
	r.makeContigStarts()
	err = r.serialize(250000, "genome.json", "data.npy")
	if err != nil {
		return err
	}
	log.Notice("Success")
	return nil
}

func (r *Anchorer) iterativeGraphMerge(paths PathSet, flanksize int64) error {
	var G Graph
	i := 0
	prevPaths := len(paths)
	graphRemake := true
	for prevPaths > 1 {
		if graphRemake {
			i++
			log.Noticef("Starting iteration %d with %d paths (L50=%d)",
				i, len(paths), getL50(paths))
			G = r.makeGraph()
		}
		CG := r.makeConfidenceGraph(G)
		paths, err := r.generatePathAndCycle(CG, flanksize)
		if err != nil {
			return err
		}
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
	return nil
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
	if smallestPath == nil {
		return nil
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
func (r *Anchorer) makeTrivialPaths(contigs []*Contig, flanksize int64) (PathSet, error) {
	// Initially make every contig a single Path object
	paths := PathSet{}
	for _, contig := range contigs {
		contig.orientation = 1
		_, err := makePath([]*Contig{contig}, paths, flanksize)
		if err != nil {
			return nil, err
		}
	}

	return paths, nil
}

// ExtractInterContigLinks extracts links from the Bamfile
func (r *Anchorer) ExtractInterContigLinks() error {
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	fh, err := os.Open(r.Bamfile)
	if err != nil {
		return err
	}
	prefix := RemoveExt(r.Bamfile)
	disfile := prefix + ".dis"
	idsfile := prefix + ".ids"

	br, _ := bam.NewReader(fh, 0)

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
		_, err := fmt.Fprintf(wids, "%s\t%d\n", ref.Name(), ref.Len())
		if err != nil {
			return err
		}
	}
	err = wids.Flush()
	if err != nil {
		return err
	}
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
		_, err := fmt.Fprintf(wdis, "%s\t%s\n", contig, arrayToString(links, ","))
		if err != nil {
			return err
		}
	}
	err = wdis.Flush()
	if err != nil {
		return err
	}
	log.Noticef("Extracted %d intra-contig and %d inter-contig links",
		intraTotal, interTotal)

	return br.Close()
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

// getLength computes the length of the path
func (r *Path) setLength() {
	r.length = 0
	for _, contig := range r.contigs {
		r.length += contig.length
	}
}

// bisect cuts the Path into two parts
// The initial implementation cuts the path into two equal halves. However, this is
// not desired for longer path since the 'internal' links should be penalized somehow.
// Therefore, we add a flanksize parameter so that we only create end nodes that are
// min(pathlength / 2, flanksize) so that we handle long paths properly
func (r *Path) bisect(flanksize int64) error {
	var contig *Contig
	var contigpos int64

	r.setLength()
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
			{0, contig.length, LNode},
		}
		contigStart += contig.length
	}
	if contig == nil {
		return fmt.Errorf("failed to find contig on the left side")
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
			{0, contig.length, RNode},
		}
	}
	return nil
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
	counts := make([]int, 0)
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
func (r *Anchorer) splitPath(piler Piler, countCutoff int, flanksize int64) (PathSet, error) {
	contigs := make([]*Contig, 0)
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
				path, err := makePath(contigs, paths, flanksize)
				if err != nil {
					return nil, err
				}
				fmt.Println(path, len(path.contigs), path.length)
				contigs = []*Contig{}
			}
		}
		contigs = append(contigs, contig)
		// fmt.Println(contig.name, contig.start, contig.orientation, strength)
	}
	// Last piece
	_, err := makePath(contigs, paths, flanksize)
	if err != nil {
		return nil, err
	}
	log.Noticef("Split into %d paths", len(paths))
	return paths, nil
}

// makePath creates a Path from contigs and set everything properly
func makePath(contigs []*Contig, paths PathSet, flanksize int64) (*Path, error) {
	path := &Path{
		contigs: contigs,
	}
	for _, contig := range contigs {
		contig.path = path
	}
	err := path.bisect(flanksize)
	if err != nil {
		return nil, err
	}
	paths[path] = true

	return path, nil
}

// findBin returns the i-th bin along the path
func findBin(contig *Contig, pos, resolution int64) int {
	position := findPosition(contig, pos)
	return int(position / resolution)
}

// serialize outputs the current path to disk
// This contains the data for jcvi.assembly.hic.heatmap()
func (r *Anchorer) serialize(res int64, jsonfile, npyfile string) error {
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
	jw := bufio.NewWriter(f)
	_, err := jw.WriteString(string(s))
	if err != nil {
		return err
	}
	err = jw.Flush()
	if err != nil {
		return err
	}
	log.Noticef("Contig stats (N=%d Length=%d) written to `%s`",
		m, r.path.length, jsonfile)

	// Serialize the pixelated matrix to NPY file
	w, _ := gonpy.NewFileWriter(npyfile)
	w.Shape = []int{m, m}
	err = w.WriteInt32(C)
	if err != nil {
		return err
	}
	log.Noticef("Matrix (resolution=%d) written to `%s`", res, npyfile)
	return f.Close()
}

// getL50 computes the L50 of all component contigs within a path
func getL50(paths PathSet) int64 {
	pathLengths := make([]int64, 0)
	for path := range paths {
		pathLengths = append(pathLengths, path.length)
	}

	return L50(pathLengths)
}

// parseTourFile parses tour file
// Only the last line is retained anc onverted into a Tour
func (r *Anchorer) parseTourFile(filename string) error {
	words, err := parseTourFile(filename)
	if err != nil {
		return err
	}
	tigs := make([]*Contig, 0)

	for _, word := range words {
		tigName, tigOrientation := word[:len(word)-1], word[len(word)-1]
		tig, ok := r.nameToContig[tigName]
		if !ok {
			log.Errorf("Contig %s not found! Skipped", tigName)
			continue
		}
		tigs = append(tigs, tig)
		tig.orientation = 1
		if tigOrientation == '-' {
			tig.orientation = -1
		}
	}

	r.path = &Path{contigs: tigs}
	r.path.setLength()
	return nil
}

// printTour logs the current tour to file
func (r *Anchorer) printTour(fwtour *os.File, label string) error {
	_, err := fwtour.WriteString(">" + label + "\n")
	if err != nil {
		return err
	}
	atoms := make([]string, len(r.path.contigs))
	for i, contig := range r.path.contigs {
		sign := "+"
		if contig.orientation < 0 {
			sign = "-"
		}
		atoms[i] = contig.name + sign
	}
	_, err = fwtour.WriteString(strings.Join(atoms, " ") + "\n")
	return err
}

// ************** Graph-related ********************

// Node is the scaffold ends, Left or Right (5` or 3`)
type Node struct {
	path   *Path // List of contigs
	sister *Node // Node of the other end
	length int64 // Typically min(pathlength / 2, flanksize)
}

// Edge is between two nodes in a graph
type Edge struct {
	a, b   *Node
	weight int64
}

// Graph is an adjacency list
type Graph map[*Node]map[*Node]int64

// nodeCmp is used often to present nodes ordered (since pointers and maps are not ordered)
func nodeCmp(a, b *Node) bool {
	return a.path.contigs[0].name < b.path.contigs[0].name ||
		(a.sister == b && a.isLNode())
}

// isLNode returns if a Node is an LNode (5`-end`)
func (r *Node) isLNode() bool {
	return r == r.path.LNode
}

// isRNode returns if a Node is an RNode (3`-end`)
func (r *Node) isRNode() bool {
	return r == r.path.RNode
}

// isReverse returns the orientation of an edge
func (r *Edge) isReverse() bool {
	return r.a.isRNode()
}

// isSister returns if the edge is internal to a contig
func (r *Edge) isSister() bool {
	return r.weight == 0
}

// makeGraph makes a contig linkage graph
func (r *Anchorer) makeGraph() Graph {
	G := Graph{}
	nIntra := 0    // becomes an intra-path link
	nInternal := 0 // internal to another path, too far away from the edge
	nUsed := 0
	// Go through the links for each node and compile edges
	for _, contig := range r.contigs {
		if contig.path == nil {
			continue
		}
		for _, link := range contig.links {
			a, b := r.linkToNodes(link)
			if a == nil || b == nil {
				nInternal++
				continue
			}
			if a == b || a.sister == b { // These links have now become intra, discard
				nIntra++
				continue
			}
			nUsed++
			r.insertEdge(G, a, b)
			r.insertEdge(G, b, a)
		}
	}

	// Normalize against the product of lengths of two paths
	for a, nb := range G {
		for b, score := range nb {
			G[a][b] = score * BigNorm / (a.length * b.length)
		}
	}

	// Print graph stats
	nEdges := 0
	for _, node := range G {
		nEdges += len(node)
	}
	nEdges /= 2 // since each edge counted twice
	log.Noticef("Graph contains %d nodes and %d edges (%d used, %d intra, %d internal)",
		len(G), nEdges, nUsed, nIntra, nInternal)
	return G
}

// makeConfidenceGraph re-calibrates the edge weight
// Steps are:
// 1 - calculate the link density as links divided by the product of two contigs
// 2 - calculate the confidence as the weight divided by the second largest edge
func (r *Anchorer) makeConfidenceGraph(G Graph) Graph {
	twoLargest := map[*Node][]int64{}

	for a, nb := range G {
		first, second := int64(0), int64(0)
		for _, score := range nb {
			if score > first {
				first, second = score, first
			} else if score < first && score > second {
				second = score
			}
		}
		twoLargest[a] = []int64{first, second}
	}
	// fmt.Println(G)

	confidenceGraph := Graph{}
	// Now a second pass to compute confidence
	for a, nb := range G {
		for b, weight := range nb {
			secondLargest := getSecondLargest(twoLargest[a], twoLargest[b])
			if secondLargest == 0 {
				continue
			}
			confidence := weight * BigNorm / secondLargest
			if confidence > BigNorm {
				if _, ok := confidenceGraph[a]; ok {
					confidenceGraph[a][b] = confidence
				} else {
					confidenceGraph[a] = map[*Node]int64{b: confidence}
				}
			}
		}
	}
	// fmt.Println(confidenceGraph)
	return confidenceGraph
}

// Get the second largest number without sorting
// a, b are both 2-item arrays
func getSecondLargest(a, b []int64) int64 {
	A := append(a, b...)
	sort.Slice(A, func(i, j int) bool {
		return A[i] < A[j]
	})
	// Some edge will appear twice in this list so need to remove it
	for i := 2; i >= 0; i-- {
		if A[i] < A[3] {
			return A[i]
		}
	}
	return A[0]
}

// getUniquePaths returns all the paths that are curerntly active
func (r *Anchorer) getUniquePaths() PathSet {
	paths := map[*Path]bool{}
	nSingletonContigs := 0
	nComplexContigs := 0
	nSingleton := 0
	nComplex := 0
	for _, contig := range r.contigs {
		path := contig.path
		if path == nil {
			continue
		}
		if len(path.contigs) == 1 {
			nSingletonContigs++
		} else {
			nComplexContigs++
		}
		if _, ok := paths[path]; ok {
			continue
		}
		paths[path] = true
		if len(path.contigs) == 1 {
			nSingleton++
		} else {
			nComplex++
		}
	}

	log.Noticef("%d paths (nComplex=%d nSingle=%d), %d contigs (nComplex=%d nSingle=%d)",
		nComplex+nSingleton, nComplex, nSingleton,
		nComplexContigs+nSingletonContigs, nComplexContigs, nSingletonContigs)

	return paths
}

// generatePathAndCycle makes new paths by merging the unique extensions
// in the graph. This first extends upstream (including the sister edge)
// and then walk downstream until it hits something seen before
func (r *Anchorer) generatePathAndCycle(G Graph, flanksize int64) (PathSet, error) {
	visited := map[*Node]bool{}
	var isCycle bool
	// We can just iterate the dictionary, however, that will not be ordered
	// we want the visit order to be stable
	orderedNodes := make([]*Node, 0)
	for a := range G {
		orderedNodes = append(orderedNodes, a)
	}
	sort.Slice(orderedNodes, func(i, j int) bool {
		return nodeCmp(orderedNodes[i], orderedNodes[j])
	})

	// Now go through all nodes
	for _, a := range orderedNodes {
		if _, ok := visited[a]; ok {
			continue
		}
		path1, path2 := make([]Edge, 0), make([]Edge, 0)
		path1, isCycle = dfs(G, a, path1, visited, true)

		if isCycle {
			path1 = breakCycle(path1)
		} else { // upstream search returns a path, we'll stitch
			delete(visited, a)
			path2, _ = dfs(G, a, path2, visited, false)
			path1 = append(reversePath(path1), path2...)
		}
		// fmt.Println("path1", path1)
		path, err := mergePath(path1, flanksize)
		if err != nil {
			return nil, err
		}
		// fmt.Println("path from", a, path)
		for _, contig := range path.contigs {
			contig.path = path
		}
	}
	return r.getUniquePaths(), nil
}

// mergePath converts a single edge path into a node path
func mergePath(path []Edge, flanksize int64) (*Path, error) {
	s := &Path{}
	for _, edge := range path {
		if !edge.isSister() {
			continue
		}
		ep := edge.a.path
		if edge.isReverse() {
			ep.reverse()
		}
		s.contigs = append(s.contigs, ep.contigs...)
	}
	err := s.bisect(flanksize)
	return s, err
}

// reversePath reverses a single edge path into its reverse direction
func reversePath(path []Edge) []Edge {
	ans := make([]Edge, 0)
	for i := len(path) - 1; i >= 0; i-- {
		ans = append(ans, Edge{
			path[i].b, path[i].a, path[i].weight,
		})
	}
	return ans
}

// breakCycle breaks a single edge path into two edge paths
// breakage occurs at the weakest link
func breakCycle(path []Edge) []Edge {
	minI, minWeight := 0, int64(math.MaxInt64)
	for i, edge := range path {
		if edge.weight > 1 && edge.weight < minWeight {
			minI, minWeight = i, edge.weight
		}
	}
	return append(path[minI+1:], path[:minI]...)
}

// dfs visits the nodes in DFS order
// Return the path and if the path is a cycle
func dfs(G Graph, a *Node, path []Edge, visited map[*Node]bool, visitSister bool) ([]Edge, bool) {
	if _, ok := visited[a]; ok { // A cycle
		return path, true
	}

	visited[a] = true
	// Alternating between sister and non-sister edges
	if visitSister {
		path = append(path, Edge{
			a, a.sister, 0,
		})
		return dfs(G, a.sister, path, visited, false)
	}
	if nb, ok := G[a]; ok {
		var maxb *Node
		maxWeight := int64(0) // for tie breaking
		for b, weight := range nb {
			if weight > maxWeight || (weight == maxWeight && nodeCmp(b, maxb)) {
				maxb, maxWeight = b, weight
			}
		}
		path = append(path, Edge{
			a, maxb, maxWeight,
		})
		return dfs(G, maxb, path, visited, true)
	}
	return path, false
}
