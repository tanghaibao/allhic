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
}

// Contig stores the name and length of each contig
type Contig struct {
	name   string
	length int
	links  []Link
}

// Link contains a specific inter-contig link
type Link struct {
	a, b       *Contig // Contig ids
	apos, bpos int     // Positions
}

// Run kicks off the merging algorithm
func (r *Anchorer) Run() {
	r.ExtractInterContigLinks()
	r.MakeGraph()
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
	total := 0
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
			continue
		}

		//         read1                                               read2
		//     ---a-- X|----- dist = a2 ----|         |--- dist = b ---|X ------ b2 ------
		//     ==============================         ====================================
		//             C1 (length L1)       |----D----|         C2 (length L2)
		// An inter-contig link
		a.links = append(a.links, Link{
			a: a, b: b, apos: apos, bpos: bpos,
		})
	}

	for _, contig := range r.contigs {
		sort.Slice(contig.links, func(i, j int) bool {
			return contig.links[i].apos < contig.links[j].apos
		})
	}

	// for i, links := range r.links {
	// 	fmt.Println(i, ":", links)
	// }

	// Write intra-links to .dis file
	for contig, links := range intraLinks {
		links = unique(links)
		total += len(links)
		fmt.Fprintf(wdis, "%s\t%s\n", contig, arrayToString(links, ","))
	}
	wdis.Flush()
	log.Noticef("Extracted %d intra-contig link groups to `%s` (total = %d)",
		len(intraLinks), disfile, total)
}

// Path is a collection of ordered contigs
type Path struct {
	contigs      []*Contig // List of contigs
	orientations []int     // 0 => +, 1 => -
	length       int       // Cumulative length of all contigs
}

// Node is the scaffold ends, Left or Right (5` or 3`)
type Node struct {
	path *Path // List of contigs
	end  int   // 0 => L, 1 => R
}

// Range tracks contig:start-end
type Range struct {
	contig *Contig
	start  int
	end    int
	node   *Node
}

// Graph is an adjacency list
type Graph map[*Node]map[*Node]float64

// Registry contains mapping from contig ID to node ID
// We iterate through 1 or 2 ranges per contig ID
type Registry map[*Contig][]Range

// LiftOver takes as input contig ID and position, returns the nodeID
func (r *Anchorer) LiftOver() {

}

// MakeGraph makes a contig linkage graph
func (r *Anchorer) MakeGraph() {
	// Initially make every contig a single Path object
	paths := make([]Path, len(r.contigs))
	nodes := make([]Node, 2*len(r.contigs))
	registry := make(Registry)
	nEdges := 0
	for i, contig := range r.contigs {
		path := Path{
			contigs:      []*Contig{contig},
			orientations: []int{0},
		}
		paths[i] = path
		path.bisect(registry, &nodes[2*i], &nodes[2*i+1])
	}
	// Go through the links for each node and compile edges
	log.Noticef("Graph contains %d nodes and %d edges", len(nodes), nEdges)
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

	// Update the registry that is used for liftOver
	for k := 0; k < i; k++ {
		contig = r.contigs[k]
		registry[contig] = []Range{
			Range{contig, 0, contig.length, LNode},
		}
	}
	registry[contig] = []Range{
		Range{contig, 0, contigpos, LNode},
		Range{contig, contigpos, contig.length, RNode},
	}
	for k := i + 1; k < len(r.contigs); k++ {
		contig = r.contigs[k]
		registry[contig] = []Range{
			Range{contig, 0, contig.length, RNode},
		}
	}
}
