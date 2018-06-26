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
	"math"
	"sort"
)

// Edge is between two nodes in a graph
type Edge struct {
	a, b   *Node
	weight float64
}

// isReverse returns the orientation of an edge
func (r *Edge) isReverse() bool {
	return r.a.end == 1
}

// isSister returns if the edge is internal to a contig
func (r *Edge) isSister() bool {
	return r.weight == 0
}

// makeGraph makes a contig linkage graph
func (r *Anchorer) makeGraph(paths []Path) Graph {
	G := make(Graph)
	r.registerPaths(paths)
	// Go through the links for each node and compile edges
	for _, contig := range r.contigs {
		for _, link := range contig.links {
			a, b := r.linkToNodes(link)
			r.insertEdge(G, a, b)
			r.insertEdge(G, b, a)
			// 	if (link.a.name == "idcChr1.ctg433" && link.b.name == "idcChr1.ctg434") ||
			// 		(link.a.name == "idcChr1.ctg434" && link.b.name == "idcChr1.ctg433") {
			// 		fmt.Println("***", link.a.name, link.b.name, link.apos, link.bpos,
			// 			a, b, G[a][b], G[b][a])
			// 	}
		}
	}
	nEdges := 0
	for _, node := range G {
		nEdges += len(node)
	}
	log.Noticef("Graph contains %d nodes and %d edges", len(G), nEdges)
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
			// TODO: Limit > 1 links here
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
func (r *Anchorer) generatePathAndCycle(G Graph) []Path {
	visited := map[*Node]bool{}
	var isCycle bool
	nodeToPath := make(map[*Node]*Path)
	paths := []Path{}
	for a := range G {
		if _, ok := visited[a]; ok {
			continue
		}
		path1, path2 := []Edge{}, []Edge{}
		path1, isCycle = dfs(G, a, path1, visited, true)
		if isCycle {
			path1 = breakCycle(path1)
			paths = append(paths, mergePath(path1, nodeToPath))
			continue
		}
		delete(visited, a)
		path2, _ = dfs(G, a, path2, visited, false)

		path1 = append(reversePath(path1), path2...)
		paths = append(paths, mergePath(path1, nodeToPath))
	}
	log.Noticef("A total of %d nodes mapped to %d paths", len(nodeToPath), len(paths))
	return paths
}

// mergePath converts a single edge path into a node path
func mergePath(path []Edge, nodeToPath map[*Node]*Path) Path {
	p := []string{}
	s := Path{}
	for _, edge := range path {
		if !edge.isSister() {
			continue
		}
		ep := edge.a.path
		tag := ""
		if edge.isReverse() {
			tag = "-"
			ep.reverse()
		}
		// TODO: take orientations into account
		s.contigs = append(s.contigs, ep.contigs...)
		s.orientations = append(s.orientations, ep.orientations...)
		nodeToPath[edge.a] = &s
		nodeToPath[edge.b] = &s

		// Special care needed for reverse orientation
		for _, contig := range ep.contigs {
			p = append(p, tag+contig.name)
		}
	}
	s.Length()
	// fmt.Println(path)
	// fmt.Println(s)
	// fmt.Println(p)
	return s
}

// reversePath reverses a single edge path into its reverse direction
func reversePath(path []Edge) []Edge {
	ans := []Edge{}
	for i := len(path) - 1; i >= 0; i-- {
		ans = append(ans, Edge{
			path[i].b, path[i].a, path[i].weight,
		})
	}
	return ans
}

// breakCycle breaks a single edge path into two edge paths
func breakCycle(path []Edge) []Edge {
	minI, minWeight := 0, math.MaxFloat64
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
