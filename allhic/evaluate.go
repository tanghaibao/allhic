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
	"math/rand"

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
	clone.Tigs = make([]Tig, r.Len())
	copy(clone.Tigs, r.Tigs)
	clone.M = r.M
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
	// fmt.Println(r.Tigs, mid)

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
			// We are looking for maximum
			score -= float64(nlinks) / dist
		}
	}
	return
}

// Mutate a Tour by applying by permutation mutation and/or splice mutation
// TODO: CHANGE TO GENOME_MUTATION (INVERSION + TRANSPOSITION)
func (r Tour) Mutate(rng *rand.Rand) {
	if rng.Float64() < 0.35 {
		gago.MutPermute(r, 3, rng)
	}
	if rng.Float64() < 0.45 {
		gago.MutSplice(r, rng)
	}
}

// Crossover a Tour with another Tour by using Partially Mixed Crossover (PMX).
func (r Tour) Crossover(q gago.Genome, rng *rand.Rand) {
	gago.CrossPMX(r, q.(Tour), rng)
}

// Clone a Tour
func (r Tour) Clone() gago.Genome {
	var clone Tour
	clone.Tigs = make([]Tig, r.Len())
	copy(clone.Tigs, r.Tigs)
	clone.M = r.M
	return clone
}

// Shuffle randomly shuffles an integer array using Knuth or Fisher-Yates
func (r Tour) Shuffle() {
	N := r.Len()
	for i := 0; i < N; i++ {
		// choose index uniformly in [i, N-1]
		j := i + rand.Intn(N-i)
		r.Tigs[j], r.Tigs[i] = r.Tigs[i], r.Tigs[j]
	}
}

// GARun set up the Genetic Algorithm and run it
func GARun(tour Tour) {
	MakeTour := func(rng *rand.Rand) gago.Genome {
		tour.Shuffle()
		c := tour.Clone()
		return c
	}

	ga := gago.Generational(MakeTour)
	ga.Initialize()
	// fmt.Println(ga.Populations)

	log.Notice("GA initialized")

	for i := 0; i < 5000; i++ {
		ga.Evolve()
		fmt.Printf("*** Generation %d ***\n", i)
		fmt.Println(ga.HallOfFame[0].Genome.(Tour).Tigs)
		fmt.Println(ga.HallOfFame[0].Fitness)
	}
}
