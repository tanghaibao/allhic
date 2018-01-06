/**
 * Filename: /Users/htang/code/allhic/allhic/orientation.go
 * Path: /Users/htang/code/allhic/allhic
 * Created Date: Friday, January 5th 2018, 4:35:57 pm
 * Author: bao
 *
 * Copyright (c) Haibao Tang
 */

package allhic

import (
	"fmt"

	"github.com/gonum/matrix/mat64"
)

// flipLog briefly logs if the orientation flip is informative
func flipLog(method string, score, scoreFlipped float64, tag string) {
	log.Noticef("%v: %.5f => %.5f %v", method, score, scoreFlipped, tag)
}

// flipAll initializes the orientations based on pairwise O matrix.
func (r *CLMFile) flipAll() {
	var (
		M mat64.Dense
		e mat64.EigenSym
	)

	N := len(r.Tigs)
	e.Factorize(r.O(), true)
	M.EigenvectorsSym(&e)
	v := M.ColView(N - 1) // v is the eigenvec corresponding to the largest eigenval
	// fmt.Printf("%0.2v\n\n", mat64.Formatted(v))

	signs := make([]byte, N)
	for i := 0; i < N; i++ {
		if v.At(i, 0) < 0 {
			signs[i] = '-'
		} else {
			signs[i] = '+'
		}
	}
	r.Signs = signs
	log.Notice("Eigenvector calculated on pairwise orientation matrix")
}

// flipWhole test flipping all contigs at the same time to see if score improves
func (r *CLMFile) flipWhole() {
	Q := r.Q()
	score := r.EvaluateQ(r.Tour, Q)
	fmt.Println(score)
}

// O yields a pairwise orientation matrix, where each cell contains the strandedness
// times the number of links between i-th and j-th contig
func (r *CLMFile) O() *mat64.SymDense {
	N := len(r.Tigs)
	P := mat64.NewSymDense(N, nil)
	for pair, contact := range r.contacts {
		score := float64(contact.strandedness * contact.nlinks)
		P.SetSym(pair.ai, pair.bi, score)
	}
	return P
}

// Q yields a contact frequency matrix when contigs are already oriented. This is a
// similar matrix as M, but rather than having the number of links in the
// cell, it points to an array that has the actual distances.
func (r *CLMFile) Q() [][]GArray {
	N := len(r.Tigs)
	P := Make2DGArraySlice(N, N)
	for pair, gdists := range r.orientedContacts {
		ai := pair.ai
		bi := pair.bi
		if r.Signs[ai] == pair.ao && r.Signs[bi] == pair.bo {
			P[ai][bi] = gdists
		}
	}
	return P
}

// EvaluateQ sums up all distance is defined as the sizes of interleaving contigs
// plus the actual link distances. Maximize Sum(1 / distance) for all links.
// For performance consideration, we actually use a histogram to approximate
// all link distances. See goldenArray() for details.
func (r *CLMFile) EvaluateQ(tour Tour, Q [][]GArray) (score float64) {
	size := tour.Len()
	cumsize := make([]int, size)
	cumSum := 0
	for i, t := range r.Tour.Tigs {
		cumsize[i] = cumSum
		cumSum += t.Size
	}
	// fmt.Println(tour.Tigs, cumsize)

	// Now add up all the pairwise scores
	for i := 0; i < size; i++ {
		a := tour.Tigs[i].Idx
		for j := i + 1; j < size; j++ {
			b := tour.Tigs[j].Idx
			dist := cumsize[j-1] - cumsize[i]
			if dist > LIMIT {
				break
			}
			for k := 0; k < BB; k++ {
				score += float64(Q[a][b][k]) / float64(GR[k]+dist)
			}
			// fmt.Println(r.Tigs[a], r.Tigs[b], Q[a][b], score)
		}
	}
	return
}
