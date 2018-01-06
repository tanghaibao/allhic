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
	fmt.Printf("%0.4v\n\n", e.Values(nil))
	fmt.Printf("%0.2v\n\n", mat64.Formatted(M.ColView(N-1)))
}

// O yields a pairwise orientation matrix, where each cell contains the strandedness
// times the number of links between i-th and j-th contig
func (r *CLMFile) O() *mat64.SymDense {
	N := len(r.Tigs)
	P := mat64.NewSymDense(N, nil)
	for _, contact := range r.contacts {
		score := float64(contact.strandedness * contact.nlinks)
		P.SetSym(contact.ai, contact.bi, score)
	}
	return P
}
