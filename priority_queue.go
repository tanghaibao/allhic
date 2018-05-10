/*
 * Filename: /Users/bao/code/allhic/priority_queue.go
 * Path: /Users/bao/code/allhic
 * Created Date: Thursday, May 10th 2018, 2:03:09 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

// A PriorityQueue implements heap.Interface and holds Items.
type PriorityQueue []*Item

// Len returns the number of items in the queu
func (pq PriorityQueue) Len() int { return len(pq) }

// Less defines the way items get ordered
func (pq PriorityQueue) Less(i, j int) bool {
	// We want Pop to give us the highest, not lowest, priority so we use greater than here.
	return pq[i].priority > pq[j].priority
}

// Swap exchanges values of two elements
func (pq PriorityQueue) Swap(i, j int) {
	pq[i], pq[j] = pq[j], pq[i]
	pq[i].index = i
	pq[j].index = j
}

// Push adds an element to the queue
func (pq *PriorityQueue) Push(x interface{}) {
	n := len(*pq)
	item := x.(*Item)
	item.index = n
	*pq = append(*pq, item)
}

// Pop removes an element with the lowest priority
func (pq *PriorityQueue) Pop() interface{} {
	old := *pq
	n := len(old)
	item := old[n-1]
	item.index = -1 // for safety
	*pq = old[0 : n-1]
	return item
}
