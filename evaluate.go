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
	"math"
	"math/rand"
	"os"
	"sort"

	"github.com/tanghaibao/gago"
)

// LIMIT determines the largest distance for two tigs to add to total score
const LIMIT = 10000000

// LimitLog is the Log of LIMIT
var LimitLog = math.Log(LIMIT)

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

// // Evaluate calculates a score for the current tour
// func (r Tour) Evaluate() (score float64) {
// 	size := r.Len()
// 	mid := make([]float64, size)
// 	cumSum := 0.0
// 	for i, t := range r.Tigs {
// 		tsize := float64(t.Size)
// 		mid[i] = cumSum + tsize/2
// 		cumSum += tsize
// 	}
// 	// Now add up all the pairwise scores
// 	for i := 0; i < size; i++ {
// 		a := r.Tigs[i].Idx
// 		for j := i + 1; j < size; j++ {
// 			b := r.Tigs[j].Idx
// 			nlinks := r.M[a][b]
// 			dist := mid[j] - mid[i]
// 			if dist > LIMIT {
// 				break
// 			}
// 			// We are looking for maximum
// 			score -= float64(nlinks) / dist
// 		}
// 	}
// 	return
// }

// Evaluate calculates a score for the current tour
func (r Tour) Evaluate() float64 {
	size := r.Len()
	mid := make([]float64, size)
	cumSum := 0.0
	for i, t := range r.Tigs {
		tsize := float64(t.Size)
		mid[i] = cumSum + tsize/2
		cumSum += tsize
	}

	score := 0.0
	// Now add up all the pairwise scores
	for i := 0; i < size; i++ {
		a := r.Tigs[i].Idx
		for j := i + 1; j < size; j++ {
			b := r.Tigs[j].Idx
			nlinks := r.M[a][b]
			dist := mid[j] - mid[i]
			if dist > LIMIT {
				break
			}
			// gago only looks at minimum =>
			// everytime we have a small dist, we reduce the total score
			// we are looking at the largest reductions from all links
			score += float64(nlinks) * (math.Log(dist) - LimitLog)
		}
	}
	return score
}

// Sample k unique integers in range [min, max) using reservoir sampling,
// specifically Vitter's Algorithm R. From gago.
func randomInts(k, min, max int, rng *rand.Rand) (ints []int) {
	ints = make([]int, k)
	for i := 0; i < k; i++ {
		ints[i] = i + min
	}
	for i := k; i < max-min; i++ {
		var j = rng.Intn(i + 1)
		if j < k {
			ints[j] = i + min
		}
	}
	return
}

// MutInversion applies inversion operation on the genome
func MutInversion(genome gago.Slice, n int, rng *rand.Rand) {
	// log.Debugf("Before MutInversion: %v", genome)
	for k := 0; k < n; k++ {
		// Choose two points on the genome
		var points = randomInts(2, 0, genome.Len(), rng)
		sort.Ints(points)
		p := points[0]
		q := points[1]
		if p == q {
			return
		}
		// Swap within range
		for i, j := p, q; i < j; i, j = i+1, j-1 {
			genome.Swap(i, j)
		}
	}
	// log.Debugf("After MutInversion: %v", genome)
}

// MutInsertion applies insertion operation on the genome
func MutInsertion(genome gago.Slice, n int, rng *rand.Rand) {
	// log.Debugf("Before MutInsertion: %v", genome)
	for k := 0; k < n; k++ {
		// Choose two points on the genome
		var points = randomInts(2, 0, genome.Len(), rng)
		p := points[0]
		q := points[1]
		if p == q {
			return
		}
		cq := genome.At(q) // Pop q and insert to p position
		if p < q {
			// Move cq to the front and push everyone right
			for i := q; i > p; i-- {
				genome.Set(i, genome.At(i-1))
			}
		} else { // q < p
			// Move cq to the back and push everyone left
			for i := q; i < p; i++ {
				genome.Set(i, genome.At(i+1))
			}
		}
		genome.Set(p, cq)
	}
	// log.Debugf("After MutInsertion: %v", genome)
}

// Mutate a Tour by applying by inversion or insertion
func (r Tour) Mutate(rng *rand.Rand) {
	if rng.Float64() < 0.5 {
		MutInversion(r, 1, rng)
	} else {
		MutInsertion(r, 1, rng)
	}
}

// Crossover a Tour with another Tour by using Partially Mixed Crossover (PMX).
func (r Tour) Crossover(q gago.Genome, rng *rand.Rand) {
	//gago.CrossPMX(r, q.(Tour), rng)
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
func (r *CLM) GARun(fwtour *os.File, npop, ngen int, mutRate, crossRate float64, phase int) Tour {
	MakeTour := func(rng *rand.Rand) gago.Genome {
		c := r.Tour.Clone()
		return c
	}

	ga := gago.GA{
		NewGenome: MakeTour,
		NPops:     1,
		PopSize:   npop,
		Model: gago.ModGenerational{
			Selector: gago.SelTournament{
				NContestants: 3,
			},
			MutRate:   mutRate,
			CrossRate: crossRate,
		},
		ParallelEval: true,
	}
	ga.Initialize()

	log.Noticef("GA initialized (npop: %v, ngen: %v, mu: %.3f, cx: %.3f)",
		npop, ngen, mutRate, crossRate)

	gen := 1
	best := -math.MaxFloat64
	updated := 0
	for ; ; gen++ {
		ga.Evolve()
		currentBest := -ga.HallOfFame[0].Fitness
		if gen%(ngen/10) == 0 {
			fmt.Printf("Current iteration GA%d-%d: max_score=%.5f\n",
				phase, gen, currentBest)
			currentBestTour := ga.HallOfFame[0].Genome.(Tour)
			r.PrintTour(fwtour, currentBestTour, fmt.Sprintf("GA%d-%d-%.5f",
				phase, gen, currentBest))
		}

		if currentBest > best {
			best = currentBest
			updated = gen
		}

		if gen-updated > ngen { // Converged
			break
		}
	}
	r.Tour = ga.HallOfFame[0].Genome.(Tour)
	return r.Tour
}
