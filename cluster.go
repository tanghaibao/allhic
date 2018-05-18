/*
 * Filename: /Users/bao/code/allhic/cluster.go
 * Path: /Users/bao/code/allhic
 * Created Date: Saturday, April 21st 2018, 4:09:48 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"container/heap"
	"fmt"
)

// Item is a generic type that stores the merges
type Item struct {
	a        int
	b        int
	priority float64
	index    int
}

// Cluster performs the hierarchical clustering
// This function is a re-implementation of the AHClustering() function in LACHESIS
func Cluster(G [][]float64, nclusters int) {

	// TODO: Skip contigs that are too small or irrelevant
	// LACHESIS also skips contigs that are thought to be centromeric

	N := len(G)
	log.Noticef("Clustering starts with %d contigs with target of %d clusters",
		N, nclusters)

	// Auxiliary data structures to facilitate cluster merging
	clusterID := make([]int, 2*N)
	clusterSize := make([]int, 2*N)
	clusterActive := make([]bool, 2*N)

	// Initially all contigs in their own cluster
	for i := 0; i < N; i++ {
		clusterID[i] = i
		clusterSize[i] = 1
		clusterActive[i] = true
	}

	// mergeScores has all possible pairwise merge scores
	var pq PriorityQueue
	for i := 0; i < N; i++ {
		for j := i + 1; j < N; j++ {
			if G[i][j] > 0 {
				pq = append(pq, &Item{
					a:        i,
					b:        j,
					priority: G[i][j],
					index:    pq.Len(),
				})
			}
		}
	}

	heap.Init(&pq)
	log.Noticef("Queue contains %d candidate merges", pq.Len())

	var item *Item
	nMerges := 0
	// The core hierarchical clustering
	for nMerges < nclusters {
		fmt.Println(pq.Len())
		// Step 1. Find the pairs of the clusters with the highest merge score
		if pq.Len() > 0 {
			item = heap.Pop(&pq).(*Item)
		}
		if pq.Len() == 0 {
			log.Notice("No more merges to do since the queue is empty")
			break
		}
		fmt.Println(item, item.priority)

		// Step 2. Merge the contig pair
		newClusterID := N + nMerges
		clusterActive[item.a] = false
		clusterActive[item.b] = false
		clusterActive[newClusterID] = true
		clusterSize[newClusterID] = clusterSize[item.a] + clusterSize[item.b]

		var newCluster []int
		for i := 0; i < N; i++ {
			if clusterID[i] == item.a || clusterID[i] == item.b {
				clusterID[i] = newClusterID
				newCluster = append(newCluster, i)
			}
		}

		// Step 3. Calculate new score entries for the new cluster
		totalLinkageByCluster := make([]float64, 2*N)
		for i := 0; i < N; i++ {
			cID := clusterID[i]
			if cID == newClusterID { // No need to calculate linkages within cluster
				continue
			}
			for _, j := range newCluster {
				totalLinkageByCluster[cID] += G[i][j]
			}
		}
		nMerges++
	}
}
