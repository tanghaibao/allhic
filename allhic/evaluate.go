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
)

// LIMIT determines the largest distance for two tigs to add to total score
const LIMIT = 10000000

// EvaluateM calculates a score for the current tour
func EvaluateM(tour Tour) (s float64) {
	size := len(tour.Tigs)
	mid := make([]float64, size)
	cumSum := 0.0
	for i, t := range tour.Tigs {
		tsize := float64(t.Size)
		mid[i] = cumSum + tsize/2
		cumSum += tsize
	}
	fmt.Println(mid)

	// Now add up all the pairwise scores
	for i := 0; i < size; i++ {
		for j := i + 1; j < size; j++ {
			a := tour.Tigs[i].Idx
			b := tour.Tigs[j].Idx
			nlinks := tour.M[a][b]
			dist := mid[j] - mid[i]
			if dist > LIMIT {
				break
			}
			s += float64(nlinks) / dist
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
