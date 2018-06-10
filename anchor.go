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
	"os"
	"sort"

	"github.com/biogo/hts/bam"
)

// Anchorer runs the merging algorithm
type Anchorer struct {
	Bamfile      string
	contigs      []*Contig
	nameToContig map[string]*Contig
	registry     Registry
}

// Contig stores the name and length of each contig
type Contig struct {
	name   string
	length int
	links  []*Link
}

// Link contains a specific inter-contig link
type Link struct {
	a, b       *Contig // Contig ids
	apos, bpos int     // Positions
}

// Run kicks off the merging algorithm
func (r *Anchorer) Run() {
	r.ExtractInterContigLinks()
	G := r.makeGraph()
	G = r.makeConfidenceGraph(G)
	r.generatePathAndCycle(G)
	log.Notice("Success")
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
	contigs      []*Contig // List of contigs
	orientations []int     // 1 => +, -1 => -
	length       int       // Cumulative length of all contigs
}

// flipOrientations reverses the orientations of all components
func (r *Path) reverse() {
	c := r.contigs
	for i, j := 0, len(c)-1; i < j; i, j = i+1, j-1 {
		c[i], c[j] = c[j], c[i]
	}
	o := r.orientations
	for i, j := 0, len(o)-1; i < j; i, j = i+1, j-1 {
		o[i], o[j] = o[j], o[i]
	}
	for i := range o {
		o[i] = -o[i]
	}
}

// Node is the scaffold ends, Left or Right (5` or 3`)
type Node struct {
	path   *Path // List of contigs
	end    int   // 0 => L, 1 => R
	sister *Node // Node of the other end
}

// Range tracks contig:start-end
type Range struct {
	start int
	end   int
	node  *Node
}

// Graph is an adjacency list
type Graph map[*Node]map[*Node]float64

// Registry contains mapping from contig ID to node ID
// We iterate through 1 or 2 ranges per contig ID
type Registry map[*Contig][]Range

// makeGraph makes a contig linkage graph
func (r *Anchorer) makeGraph() Graph {
	// Initially make every contig a single Path object
	paths := make([]Path, len(r.contigs))
	nodes := make([]Node, 2*len(r.contigs))
	r.registry = make(Registry)
	for i, contig := range r.contigs {
		path := Path{
			contigs:      []*Contig{contig},
			orientations: []int{1},
		}
		paths[i] = path
		path.bisect(r.registry, &nodes[2*i], &nodes[2*i+1])
	}

	G := make(Graph)
	// Go through the links for each node and compile edges
	for _, contig := range r.contigs {
		for _, link := range contig.links {
			a, b := r.linkToNodes(link)
			r.insertEdge(G, a, b)
			r.insertEdge(G, b, a)
		}
	}
	nEdges := 0
	for _, node := range G {
		nEdges += len(node)
	}
	log.Noticef("Graph contains %d nodes and %d edges", len(G), nEdges)
	// fmt.Println(G)
	return G
}

// contigToNode takes as input contig and position, returns the nodeID
func (r *Anchorer) contigToNode(contig *Contig, pos int) *Node {
	for _, rr := range r.registry[contig] {
		if pos >= rr.start && pos < rr.end {
			return rr.node
		}
	}
	return nil
}

// linkToNodes takes as input link, returns two nodeIDs
func (r *Anchorer) linkToNodes(link *Link) (*Node, *Node) {
	a := r.contigToNode(link.a, link.apos)
	b := r.contigToNode(link.b, link.bpos)
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

// Length returns the cumulative length of all contigs
func (r *Path) Length() int {
	r.length = 0
	for _, contig := range r.contigs {
		r.length += contig.length
	}
	return r.length
}

// findMidPoint find the center of the a path of contigs
func (r *Path) findMidPoint() (int, int) {
	midpos := r.Length() / 2
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
	if r.orientations[i] == -1 {
		contigpos = contig.length - contigpos
	}
	return i, contigpos
}

// bisect cuts the Path into two parts
func (r *Path) bisect(registry Registry, LNode, RNode *Node) {
	var contig *Contig
	i, contigpos := r.findMidPoint()
	contig = r.contigs[i]
	// L node
	*LNode = Node{
		path: r,
		end:  0,
	}
	// R node
	*RNode = Node{
		path: r,
		end:  1,
	}
	LNode.sister = RNode
	RNode.sister = LNode

	// Update the registry to convert contig:start-end range to nodes
	for k := 0; k < i; k++ { // Left contigs
		contig = r.contigs[k]
		registry[contig] = []Range{
			Range{0, contig.length, LNode},
		}
	}

	// Handles the middle contig as a special case
	var leftRange, rightRange Range
	if r.orientations[i] > 0 { // Forward orientation
		leftRange = Range{0, contigpos, LNode}
		rightRange = Range{contigpos, contig.length, RNode}
	} else { // Reverse orientation
		leftRange = Range{0, contigpos, RNode}
		rightRange = Range{contigpos, contig.length, LNode}
	}
	registry[contig] = []Range{
		leftRange, rightRange,
	}

	for k := i + 1; k < len(r.contigs); k++ { // Right contigs
		contig = r.contigs[k]
		registry[contig] = []Range{
			Range{0, contig.length, RNode},
		}
	}
}
