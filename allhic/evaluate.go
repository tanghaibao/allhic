/**
 * Filename: /Users/bao/code/allhic/allhic/evaluate.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Wednesday, January 3rd 2018, 9:40:36 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"fmt"

	"github.com/MaxHalford/gago"
)

// LIMIT determines the largest distance for two tigs to add to total score
const LIMIT = 10000000

// We will implement the Slice interface here, key ideas borrowed from:
// https://github.com/MaxHalford/gago-examples/blob/master/tsp_grid/main.go

// At method from Slice
func (r Tour) At(i int) interface{} {
	return r.Tigs[i]
}

// Set method from Slice
func (r Tour) Set(i int, v interface{}) {
	r.Tigs[i] = v.(Tig)
}

// Len method from Slice
func (r Tour) Len() int {
	return len(r.Tigs)
}

// Swap method from Slice
func (r Tour) Swap(i, j int) {
	r.Tigs[i], r.Tigs[j] = r.Tigs[j], r.Tigs[i]
}

// Slice method from Slice
func (r Tour) Slice(a, b int) gago.Slice {
	return Tour{r.Tigs[a:b], r.M}
}

// Split method from Slice
func (r Tour) Split(k int) (gago.Slice, gago.Slice) {
	return Tour{r.Tigs[:k], r.M}, Tour{r.Tigs[k:], r.M}
}

// Append method from Slice
func (r Tour) Append(q gago.Slice) gago.Slice {
	return Tour{append(r.Tigs, q.(Tour).Tigs...), r.M}
}

// Replace method from Slice
func (r Tour) Replace(q gago.Slice) {
	copy(r.Tigs, q.(Tour).Tigs)
}

// Copy method from Slice
func (r Tour) Copy() gago.Slice {
	var clone Tour
	copy(clone.Tigs, r.Tigs)
	return clone
}

// Evaluate calculates a score for the current tour
func (r Tour) Evaluate() (score float64) {
	size := r.Len()
	mid := make([]float64, size)
	cumSum := 0.0
	for i, t := range r.Tigs {
		tsize := float64(t.Size)
		mid[i] = cumSum + tsize/2
		cumSum += tsize
	}
	fmt.Println(mid)

	// Now add up all the pairwise scores
	for i := 0; i < size; i++ {
		for j := i + 1; j < size; j++ {
			a := r.Tigs[i].Idx
			b := r.Tigs[j].Idx
			nlinks := r.M[a][b]
			dist := mid[j] - mid[i]
			if dist > LIMIT {
				break
			}
			score += float64(nlinks) / dist
		}
	}
	return
}

// MakeVector returns a random vector by generating 5 values uniformally
// distributed between -10 and 10.
// func MakeVector(rng *rand.Rand) gago.Genome {
// 	return Vector(gago.InitUnifFloat64(2, -20, 20, rng))
// }

// GASetup set up the Genetic algorithm
// func GASetup() {
// 	var ga = gago.Generational(MakeVector)
// 	ga.Initialize()
// }
