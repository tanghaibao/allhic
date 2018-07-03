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
func (r *Anchorer) makeGraph(paths PathSet) Graph {
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
				// There can be ties in terms of scores
				// if len(confidenceGraph[a]) > 1 {
				// 	// if confidence < BigNorm*101/100 {
				// 	fmt.Println(a.path, "<=>", b.path, a.isLNode(), b.isLNode(),
				// 		G[a], "\n", confidenceGraph[a],
				// 		twoLargest[a], twoLargest[b], secondLargest)
				// 	for k, v := range confidenceGraph[a] {
				// 		fmt.Println(k, k.path, k.sister, v)
				// 	}
				// }
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
func (r *Anchorer) generatePathAndCycle(G Graph, flanksize int64) PathSet {
	visited := map[*Node]bool{}
	var isCycle bool
	var path *Path
	// We can just iterate the dictionary, however, that will not be ordered
	// we want the visit order to be stable
	orderedNodes := []*Node{}
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
		path1, path2 := []Edge{}, []Edge{}
		path1, isCycle = dfs(G, a, path1, visited, true)

		if isCycle {
			path1 = breakCycle(path1)
		} else { // upstream search returns a path, we'll stitch
			delete(visited, a)
			path2, _ = dfs(G, a, path2, visited, false)
			path1 = append(reversePath(path1), path2...)
		}
		// fmt.Println("path1", path1)
		path = mergePath(path1, flanksize)
		// fmt.Println("path from", a, path)
		for _, contig := range path.contigs {
			contig.path = path
		}
	}
	return r.getUniquePaths()
}

// mergePath converts a single edge path into a node path
func mergePath(path []Edge, flanksize int64) *Path {
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
	s.bisect(flanksize)
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
// breakage occurs at the weakest link
func breakCycle(path []Edge) []Edge {
	minI, minWeight := 0, int64(math.MaxInt64)
	var minEdge *Edge
	for i, edge := range path {
		if edge.weight > 1 && (edge.weight < minWeight ||
			(edge.weight == minWeight && nodeCmp(edge.a, minEdge.a))) {
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
		// if len(nb) > 1 {
		// 	fmt.Println(a, nb, b, maxWeight)
		// }
		path = append(path, Edge{
			a, maxb, maxWeight,
		})
		return dfs(G, maxb, path, visited, true)
	}
	return path, false
}
