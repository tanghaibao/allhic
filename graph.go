/*
 * Filename: /Users/bao/code/allhic/graph.go
 * Path: /Users/bao/code/allhic
 * Created Date: Monday, June 4th 2018, 11:37:27 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"fmt"
	"sort"
)

// MakeGraph makes a contig linkage graph
func (r *Anchorer) MakeGraph() Graph {
	// Initially make every contig a single Path object
	paths := make([]Path, len(r.contigs))
	nodes := make([]Node, 2*len(r.contigs))
	r.registry = make(Registry)
	for i, contig := range r.contigs {
		path := Path{
			contigs:      []*Contig{contig},
			orientations: []int{0},
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

// calculateEdges re-calibrates the edge weight
// Steps are:
// 1 - calculate the link density as links divided by the product of two contigs
// 2 - calculate the confidence as the weight divided by the second largest edge
func (r *Anchorer) calculateEdges(G Graph) {
	twoLargest := map[*Node][]float64{}
	for a, nb := range G {
		first, second := 0.0, 0.0
		for b, score := range nb {
			score /= float64(a.path.length) * float64(b.path.length)
			if score > first {
				first, second = score, first
			} else if score <= first && score > second {
				second = score
			}
			G[a][b] = score
		}
		twoLargest[a] = []float64{first, second}
	}
	// Now a second pass to compute confidence
	for a, nb := range G {
		for b := range nb {
			secondLargest := getSecondLargest(twoLargest[a], twoLargest[b])
			if G[a][b]/secondLargest > 1 {
				fmt.Println(a, b, G[a][b]/secondLargest, twoLargest[a], twoLargest[b])
			}
			G[a][b] /= secondLargest
		}
	}
}

// Get the second largest number without sorting
// a, b are both 2-item arrays
func getSecondLargest(a, b []float64) float64 {
	A := append(a, b...)
	sort.Float64s(A)
	largest, secondLargest := A[3], A[2]
	// Some edge will appear twice in this list so need to remove it
	if largest == secondLargest && A[1] > 0 { // Is precision an issue here?
		secondLargest = A[1]
	}

	return secondLargest
}
