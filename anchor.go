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
	contigs      []Contig
	nameToContig map[string]int
	links        [][]Link
}

// Contig stores the name and length of each contig
type Contig struct {
	id     int
	name   string
	length int
}

// Link contains a specific inter-contig link
type Link struct {
	a, b       int // Contig ids
	apos, bpos int // Positions
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

	r.nameToContig = make(map[string]int)
	refs := br.Header().Refs()
	for _, ref := range refs {
		contig := Contig{
			id:     len(r.contigs),
			name:   ref.Name(),
			length: ref.Len(),
		}
		r.contigs = append(r.contigs, contig)
		r.nameToContig[contig.name] = contig.id
		fmt.Fprintf(wids, "%s\t%d\n", ref.Name(), ref.Len())
	}
	wids.Flush()
	log.Noticef("Extracted %d contigs to `%s`", len(r.contigs), idsfile)

	// Import links into pairs of contigs
	total := 0
	intraLinks := make(map[string][]int)
	r.links = make([][]Link, len(r.contigs))
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

		r.links[a] = append(r.links[a], Link{
			a: a, b: b, apos: apos, bpos: bpos,
		})
	}

	for _, links := range r.links {
		sort.Slice(links, func(i, j int) bool {
			return links[i].apos < links[j].apos
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
	contigs      []int // List of contigs
	orientations []int // 0 => +, 1 => -
	length       int   // Cumulative length of all contigs
}

// Node is the scaffold ends, Left or Right (5` or 3`)
type Node struct {
	path  *Path // List of contigs
	end   int   // 0 => L, 1 => R
	links []Link
}

// Graph is an adjacency list
type Graph map[*Node]map[*Node]float64

// LiftOver takes as input contig ID and position, returns the nodeID
func (r *Anchorer) LiftOver() {

}

// MakeGraph makes a contig linkage graph
func (r *Anchorer) MakeGraph() {
	// Initially make every contig a single Path object
	paths := make([]Path, len(r.contigs))
	nodes := make([]Node, 2*len(r.contigs))
	nEdges := 0
	for i, contig := range r.contigs {
		path := Path{
			contigs:      []int{i},
			orientations: []int{0},
			length:       contig.length,
		}
		paths[i] = path
		bLinks, eLinks := bisect(contig.length/2, r.links[i])
		// B node
		nodes[2*i] = Node{
			path:  &path,
			end:   0,
			links: bLinks,
		}
		// E node
		nodes[2*i+1] = Node{
			path:  &path,
			end:   1,
			links: eLinks,
		}
	}
	// Go through the links for each node and compile edges
	log.Noticef("Graph contains %d nodes and %d edges", len(nodes), nEdges)
}

func bisect(mid int, links []Link) ([]Link, []Link) {
	// When links are sorted, we simply perform a binary search
	// to find the midpoint
	midID := sort.Search(len(links), func(i int) bool {
		return links[i].apos > mid
	})
	return links[:midID], links[midID:]
}
