/**
 * Filename: /Users/htang/code/allhic/script/newman_modularity.go
 * Path: /Users/htang/code/allhic/script
 * Created Date: Thursday, January 25th 2018, 6:01:39 pm
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

// "Prototype of the Newman modularity inference method. Implemented based on
// Newman. (2006) Modularity and community structure in networks.

package allhic

import (
	"fmt"
	"sync"

	"github.com/gonum/matrix/mat64"
)

// Edge represents an edge in the graph
type Edge struct {
	a, b   int
	weight int
}

// Graph represents the graph
type Graph struct {
	A        [][]int
	nodes    []string
	nodesIdx map[string]int
	m, n     int
}

// NewmanSubPartition further split a partition into smaller parts
// B* = B_ij - d_ij * Sum_k belongs to g B_ik
// d_ij is the Kronecker delta function: d_ij = 0 when i != j, 1 otherwise
func NewmanSubPartition(g Graph, B *mat64.SymDense, selected []int) [][]int {
	var ans [][]int

	n := len(selected)
	// B*: modularity matrix for partition G only
	Bs := mat64.NewSymDense(n, nil)
	// Map the original idx to the index within this partition
	index := make(map[int]int)
	for i, idx := range selected {
		index[idx] = i
	}
	for i := 0; i < n; i++ {
		Cii := 0.0
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			b := B.At(selected[i], selected[j])
			Bs.SetSym(i, j, b)
			Cii -= b
		}
		Bs.SetSym(i, i, Cii)
	}

	score, s := GetPartition(g.m, n, Bs)
	if score > EPS {
		log.Noticef("Final Q = %.5f", score)
	} else {
		log.Notice("Cannot further divide this graph")
		return append(ans, selected)
	}

	// Get into partitions
	var partA, partB []int
	for i, num := range s {
		if num > 0 {
			partA = append(partA, selected[i])
		} else {
			partB = append(partB, selected[i])
		}
	}
	fmt.Println(s)
	fmt.Println(partA)
	fmt.Println(partB)

	for _, part := range NewmanSubPartition(g, B, partA) {
		ans = append(ans, part)
	}
	for _, part := range NewmanSubPartition(g, B, partB) {
		ans = append(ans, part)
	}
	return ans
}

// NewmanPartition partitions the graph to maximize 'modularity'
// Newman modularity inference method:
//
//     Q = 1/4m Sum_ij(A_ij - k_i * k_j / 2m) s_i * s_j
//
//     m: total number of edges
//     A_ij: adjancency matrix
//     k_i, k_j: degree of node i, j
//     s_i, s_j: partition of node i, j, either +1 or -1
//
//     We can conveniently write Q in matrix form
//
//     Q = 1/4m s.T * B * s
//
//     Where B_ij = A_ij - k_i * k_j / 2m
func NewmanPartition(g Graph) [][]int {
	var ans [][]int

	B := mat64.NewSymDense(g.n, nil)
	k := make([]int, g.n)
	for i := 0; i < g.n; i++ {
		for j := 0; j < g.n; j++ {
			k[i] += g.A[i][j]
		}
	}
	// fmt.Println(k)

	for i := 0; i < g.n; i++ {
		for j := i; j < g.n; j++ {
			B.SetSym(i, j, float64(g.A[i][j])-float64(k[i]*k[j])/float64(2*g.m))
		}
	}

	score, s := GetPartition(g.m, g.n, B)
	if score > EPS {
		log.Noticef("Final Q = %.5f", score)
	} else {
		log.Notice("Cannot further divide this graph")
		selected := make([]int, g.n)
		for i := 0; i < g.n; i++ {
			selected[i] = i
		}
		return append(ans, selected)
	}

	// Get into partitions
	var partA, partB []int
	for i, num := range s {
		if num > 0 {
			partA = append(partA, i)
		} else {
			partB = append(partB, i)
		}
	}
	fmt.Println(s)
	fmt.Println(partA)
	fmt.Println(partB)
	for _, part := range NewmanSubPartition(g, B, partA) {
		ans = append(ans, part)
	}
	for _, part := range NewmanSubPartition(g, B, partB) {
		ans = append(ans, part)
	}
	return ans
}

// EvaluateQ calculates the Q score
func EvaluateQ(s []int, m, n int, B *mat64.SymDense) float64 {
	ans := 0.0
	for i := 0; i < n; i++ {
		ans += B.At(i, i)
		for j := i + 1; j < n; j++ {
			ans += 2 * B.At(i, j) * float64(s[i]*s[j])
		}
	}
	return ans / float64(4*m)
}

// EvaluateDeltaQ avoids recomputing the Q score
func EvaluateDeltaQ(s []int, m, n int, B *mat64.SymDense, i int) float64 {
	ans := 0.0
	// We flip the partition of s[i]
	for j := 0; j < n; j++ {
		// - s_i * s_j + (- s_i * s_j)
		if i != j {
			ans -= 4 * B.At(i, j) * float64(s[i]*s[j])
		}
	}
	return ans / float64(4*m)
}

// Score of a particular signs array
type Score struct {
	score float64
	idx   int
}

// GetPartition returns score and partition s
func GetPartition(m, n int, B *mat64.SymDense) (float64, []int) {
	var (
		M mat64.Dense
		e mat64.EigenSym
	)
	// Eigen decomposition
	e.Factorize(B, true)
	M.EigenvectorsSym(&e)
	v := M.ColView(n - 1) // Eigenvector corresponding to the largest eigenval
	// fmt.Printf("%0.2v\n\n", mat64.Formatted(v))

	s := make([]int, n)
	for i := 0; i < n; i++ {
		if v.At(i, 0) < 0 {
			s[i] = -1
		} else {
			s[i] = 1
		}
	}

	// Refinement
	return RefinePartition(s, m, n, B)
}

// RefinePartition refines partition by testing if flipping partition for each
// contig could increase score Q. At each iteration, the largest deltaQ is selected.
func RefinePartition(s []int, m, n int, B *mat64.SymDense) (float64, []int) {
	var wg sync.WaitGroup

	origScore := EvaluateQ(s, m, n, B)
	for {
		ch := make(chan Score, n)
		for i := 0; i < n; i++ {
			wg.Add(1)
			go func(i int) {
				defer wg.Done()
				newScore := EvaluateDeltaQ(s, m, n, B, i)
				ch <- Score{newScore, i}
			}(i)
		}

		// Wait for all workers to finish
		wg.Wait()
		close(ch)

		// Find the best score
		best := Score{-1.0, -1}
		for e := range ch {
			if e.score > best.score {
				best = e
			}
		}

		if best.score > 0 {
			log.Noticef("ACCEPTED: Q = %.5f, Q' = %.5f (flip %d)",
				origScore, origScore+best.score, best.idx)
			s[best.idx] = -s[best.idx]
			origScore += best.score
		} else {
			break
		}
	}

	return origScore, s
}
