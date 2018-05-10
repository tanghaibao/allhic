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
	at       int
	bt       int
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
	contigToClusterID := make([]int, 2*N)
	clusterSize := make([]int, 2*N)
	clusterActive := make([]bool, 2*N)

	// Initially all contigs in their own cluster
	for i := 0; i < N; i++ {
		contigToClusterID[i] = i
		clusterSize[i] = 1
		clusterActive[i] = true
	}

	// mergeScores has all possible pairwise merge scores
	var pq PriorityQueue
	for i := 0; i < N; i++ {
		for j := i + 1; j < N; j++ {
			if G[i][j] > 0 {
				pq = append(pq, &Item{
					at:       i,
					bt:       j,
					priority: G[i][j],
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
		}
		fmt.Println(item, item.priority)
		nMerges++
	}
}
