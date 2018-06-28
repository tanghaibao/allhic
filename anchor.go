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
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strings"

	"github.com/biogo/hts/bam"
)

// Anchorer runs the merging algorithm
type Anchorer struct {
	Bamfile      string
	contigs      []*Contig
	nameToContig map[string]*Contig
}

// Contig stores the name and length of each contig
type Contig struct {
	name        string
	length      int
	links       []*Link
	path        *Path
	orientation int // 1 => +, -1 => -
	segments    []Range
}

// Link contains a specific inter-contig link
type Link struct {
	a, b       *Contig // Contig ids
	apos, bpos int     // Positions
}

// Run kicks off the merging algorithm
func (r *Anchorer) Run() {
	var G Graph
	r.ExtractInterContigLinks()
	paths := r.makeTrivialPaths(r.contigs)
	prevPaths := len(paths)
	i := 0
	graphRemake := true
	for prevPaths > 1 {
		if graphRemake {
			i++
			log.Noticef("Starting iteration %d with %d paths", i, len(paths))
			G = r.makeGraph(paths)
		}
		CG := r.makeConfidenceGraph(G)
		paths = r.generatePathAndCycle(CG)
		// printPaths(paths)
		if len(paths) == prevPaths {
			paths = r.removeSmallestPath(paths, G)
			graphRemake = true
		}
		prevPaths = len(paths)
	}

	// Test split the final path
	res, d := 500000, 4
	r.splitPath(paths[0], res, d)
	log.Notice("Success")
}

// removeSmallestPath forces removal of the smallest path so that we can continue
// with the merging. This is also a small operation so we'll have to just modify
// the graph only slightly
func (r *Anchorer) removeSmallestPath(paths []*Path, G Graph) []*Path {
	smallestPath := paths[0]
	for _, path := range paths {
		if path.length < smallestPath.length {
			smallestPath = path
		}
	}

	// Un-assign the contigs
	for _, contig := range smallestPath.contigs {
		contig.path = nil
	}
	// // Inactivate the nodes
	// for _, node := range smallestPath.nodes {
	// 	if nb, ok := G[node]; ok {
	// 		delete(nb, node)
	// 	}
	// 	delete(G, node)
	// }

	return r.getUniquePaths()
}

// printPaths shows the current details of the clustering
func printPaths(paths []*Path) {
	for _, path := range paths {
		fmt.Println(path)
	}
}

// makeTrivialPaths starts the initial construction of Path object, with one
// contig per Path (trivial Path)
func (r *Anchorer) makeTrivialPaths(contigs []*Contig) []*Path {
	// Initially make every contig a single Path object
	paths := make([]*Path, len(contigs))
	for i, contig := range contigs {
		contig.orientation = 1
		paths[i] = &Path{
			contigs: []*Contig{contig},
		}
		paths[i].computeLength()
		contig.path = paths[i]
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
			length: ref.Len(),
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
			a: a, b: b, apos: apos, bpos: bpos,
		})
		interTotal++
	}

	for _, contig := range r.contigs {
		sort.Slice(contig.links, func(i, j int) bool {
			return contig.links[i].apos < contig.links[j].apos
		})
		// Sanity check - this produced inconsistent pairing
		// if contig.name == "idcChr1.ctg434" {
		// 	for _, link := range contig.links {
		// 		if link.b.name == "idcChr1.ctg433" {
		// 			fmt.Println(link.a.name, link.b.name, link.apos, link.bpos)
		// 		}
		// 	}
		// }
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

// Path is a collection of ordered contigs
type Path struct {
	contigs []*Contig // List of contigs
	nodes   [2]*Node  // Two nodes at each end
	length  int       // Cumulative length of all contigs
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

// Range tracks contig:start-end
type Range struct {
	start int
	end   int
	node  *Node
}

// Registry contains mapping from contig ID to node ID
// We iterate through 1 or 2 ranges per contig ID
type Registry map[*Contig][]Range

// contigToNode takes as input contig and position, returns the nodeID
func contigToNode(contig *Contig, pos int) *Node {
	for _, rr := range contig.segments { // multiple 'segments'
		if pos >= rr.start && pos < rr.end {
			return rr.node
		}
	}
	log.Errorf("%s:%d not found", contig.name, pos)
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
		G[a] = map[*Node]float64{b: 1}
	}
}

// computeLength returns the cumulative length of all contigs
func (r *Path) computeLength() int {
	r.length = 0
	for _, contig := range r.contigs {
		r.length += contig.length
	}
	return r.length
}

// findMidPoint find the center of a path for bisect
func (r *Path) findMidPoint() (int, int) {
	midpos := r.length / 2
	cumsize := 0
	i := 0
	var contig *Contig
	for i, contig = range r.contigs {
		// midpos must be within this contig
		if cumsize+contig.length > midpos {
			break
		}
		cumsize += contig.length
	}
	contigpos := midpos - cumsize

	// --> ----> <-------- ------->
	//              | mid point here
	if contig.orientation == -1 {
		contigpos = contig.length - contigpos
	}
	return i, contigpos
}

// bisect cuts the Path into two parts
func (r *Path) bisect() {
	var contig *Contig
	i, contigpos := r.findMidPoint()
	contig = r.contigs[i]

	// L node
	LNode := &Node{
		path: r,
		end:  0,
	}
	// R node
	RNode := &Node{
		path: r,
		end:  1,
	}
	r.nodes[0] = LNode
	r.nodes[1] = RNode
	LNode.sister = RNode
	RNode.sister = LNode

	// Update the registry to convert contig:start-end range to nodes
	for k := 0; k < i; k++ { // Left contigs
		contig = r.contigs[k]
		contig.segments = []Range{
			Range{0, contig.length, LNode},
		}
	}

	// Handles the middle contig as a special case
	var leftRange, rightRange Range
	contig = r.contigs[i]
	if contig.orientation > 0 { // Forward orientation
		leftRange = Range{0, contigpos, LNode}
		rightRange = Range{contigpos, contig.length, RNode}
	} else { // Reverse orientation
		leftRange = Range{0, contigpos, RNode}
		rightRange = Range{contigpos, contig.length, LNode}
	}
	contig.segments = []Range{
		leftRange, rightRange,
	}

	for k := i + 1; k < len(r.contigs); k++ { // Right contigs
		contig = r.contigs[k]
		contig.segments = []Range{
			Range{0, contig.length, RNode},
		}
	}
}

// SparseMatrix stores a big square matrix that is sparse
type SparseMatrix []map[int]int

// findBin returns the i-th bin along the path
func findBin(contigStarts map[*Contig]int, contig *Contig, pos, resolution int) int {
	start := contigStarts[contig]
	offset := pos
	if contig.orientation < 0 {
		offset = contig.length - pos
	}
	return (start + offset) / resolution
}

// splitPath takes a path and look at joins that are weak
// Scans the path at certain resolution r, and search radius is d
func (r *Anchorer) splitPath(path *Path, res, d int) {
	contigStarts := map[*Contig]int{}
	pos := 0
	for _, contig := range path.contigs {
		contigStarts[contig] = pos
		pos += contig.length
	}

	// Look at all intra-path links, then store the counts to a
	// sparse matrix, indexed by i, j, C[i, j] = # of links between
	// i-th locus and j-th locus
	bins := int(math.Ceil(float64(path.length) / float64(res)))
	log.Noticef("Contains %d bins at resolution %d bp", bins, res)
	// Initialize the sparse matrix
	C := make(SparseMatrix, bins)
	for i := 0; i < bins; i++ {
		C[i] = map[int]int{}
	}

	for _, contig := range path.contigs {
		for _, link := range contig.links {
			if _, ok := contigStarts[link.b]; !ok {
				continue
			}
			a := findBin(contigStarts, link.a, link.apos, res)
			b := findBin(contigStarts, link.b, link.bpos, res)
			if _, ok := C[a][b]; ok {
				C[a][b]++
			} else {
				C[a][b] = 1
			}

			if _, ok := C[b][a]; ok {
				C[b][a]++
			} else {
				C[b][a] = 1
			}
		}
	}
	printSparseMatrix(C, d)
}

// printMatrix shows all the entries in the matrix C that are higher than a certain
// cutoff, like 95-th percentile of all cells
func printSparseMatrix(C SparseMatrix, d int) {
	values := []int{}
	for a := range C {
		for _, val := range C[a] {
			values = append(values, val)
		}
	}

	sort.Ints(values)
	cutoff := values[len(values)*95/100]
	log.Noticef("Cutoff of cell value is at %d", cutoff)
	scores := []float64{}
	for a := range C {
		score := scoreTriangle(C, a, d, cutoff)
		scores = append(scores, score)
	}
	fmt.Println(scores)
}

// scoreTriangle sums up all the cells in the 1st quadrant that are d-distance
// away from the diagonal
func scoreTriangle(C SparseMatrix, a, d, cutoff int) float64 {
	expected := 0
	score := 0
	for i := a + 1; i < len(C); i++ {
		for j := a - 1; j >= 0 && i-j <= d; j-- {
			if _, ok := C[i][j]; ok {
				score += min(C[i][j], cutoff)
			}
			expected += cutoff
		}
	}
	fmt.Println(a, score, expected)

	// We are interested in finding all the misjoins, misjoins are dips in the
	// link coverage
	return float64(score) / float64(expected)
}
