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

// Edge is between two nodes in a graph
type Edge struct {
	a, b   *Node
	weight float64
}

// MakeGraph makes a contig linkage graph
func (r *Anchorer) makeGraph() Graph {
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
func (r *Anchorer) makeConfidenceGraph(G Graph) Graph {
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

	confidenceGraph := make(Graph)
	// Now a second pass to compute confidence
	for a, nb := range G {
		for b := range nb {
			secondLargest := getSecondLargest(twoLargest[a], twoLargest[b])
			G[a][b] /= secondLargest
			if G[a][b] > 1 {
				if _, aok := confidenceGraph[a]; aok {
					confidenceGraph[a][b] = G[a][b]
				} else {
					confidenceGraph[a] = map[*Node]float64{b: G[a][b]}
				}
			}
		}
	}
	return confidenceGraph
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

// generatePathAndCycle makes new paths by merging the unique extensions
// in the graph
func (r *Anchorer) generatePathAndCycle(G Graph) {
	fmt.Println(G)
	visited := map[*Node]bool{}
	var isCycle bool
	// paths := [][]*Node{}
	for a := range G {
		if _, ok := visited[a]; ok {
			continue
		}
		sisterPath := []Edge{}
		sisterPath, isCycle = dfs(G, a, sisterPath, visited, true)
		if isCycle {
			fmt.Println(sisterPath)
			continue
		}
		delete(visited, a)
		nonSisterPath := []Edge{}
		nonSisterPath, _ = dfs(G, a, nonSisterPath, visited, false)
		// fmt.Println(sisterPath)
		// fmt.Println(nonSisterPath)
	}
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
		var b *Node
		var weight float64
		for b, weight = range nb {
			break
		}
		path = append(path, Edge{
			a, b, weight,
		})
		return dfs(G, b, path, visited, true)
	}
	return path, false
}
