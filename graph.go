/*
 * Filename: /Users/bao/code/allhic/graph.go
 * Path: /Users/bao/code/allhic
 * Created Date: Monday, June 4th 2018, 11:37:27 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

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
	for a, nb := range G {
		for b := range nb {
			G[a][b] /= float64(a.path.length) * float64(b.path.length)
		}
	}
	// secondLargest := make(map[*Node]float64)
}
