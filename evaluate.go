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

	"github.com/MaxHalford/eaopt"
)

// LIMIT determines the largest distance for two tigs to add to total score
const LIMIT = 10000000

// LimitLog is the Log of LIMIT
var LimitLog = math.Log(LIMIT)

// We will implement the Slice interface here, key ideas borrowed from:
// https://github.com/MaxHalford/eaopt-examples/blob/master/tsp_grid/main.go

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
func (r Tour) Slice(a, b int) eaopt.Slice {
	return Tour{r.Tigs[a:b], r.M}
}

// Split method from Slice
func (r Tour) Split(k int) (eaopt.Slice, eaopt.Slice) {
	return Tour{r.Tigs[:k], r.M}, Tour{r.Tigs[k:], r.M}
}

// Append method from Slice
func (r Tour) Append(q eaopt.Slice) eaopt.Slice {
	return Tour{append(r.Tigs, q.(Tour).Tigs...), r.M}
}

// Replace method from Slice
func (r Tour) Replace(q eaopt.Slice) {
	copy(r.Tigs, q.(Tour).Tigs)
}

// Copy method from Slice
func (r Tour) Copy() eaopt.Slice {
	var clone Tour
	clone.Tigs = make([]Tig, r.Len())
	copy(clone.Tigs, r.Tigs)
	clone.M = r.M
	return clone
}

// EvaluateSumLog calculates a score for the current tour
func (r Tour) EvaluateSumLog() (float64, error) {
	//func (r Tour) Evaluate() (float64, error) {
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
			// This serves two purposes:
			// 1. Break earlier reduces the amount of calculation
			// 2. Ignore distant links so that telomeric regions don't come
			//    to be adjacent (based on Ler0 data)
			if dist > LIMIT {
				break
			}
			// eaopt only looks at minimum =>
			// everytime we have a small dist, we reduce the total score
			// we are looking at the largest reductions from all links
			score += float64(nlinks) * (math.Log(dist) - LimitLog)
		}
	}
	return score, nil
}

// Evaluate calculates a score for the current tour
func (r Tour) Evaluate() (float64, error) {
	//func (r Tour) EvaluateSumRecip() (float64, error) {
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
			// We are looking for maximum
			score -= float64(nlinks) / dist
		}
	}
	return score, nil
}

// randomTwoInts is a faster version than randomInts above
func randomTwoInts(genome eaopt.Slice, rng *rand.Rand) (int, int) {
	n := genome.Len()
	p := rng.Intn(n)
	q := rng.Intn(n)
	if p > q {
		p, q = q, p
	}
	return p, q
}

// MutInversion applies inversion operation on the genome
func MutInversion(genome eaopt.Slice, rng *rand.Rand) {
	// log.Debugf("Before MutInversion: %v", genome)
	// Choose two points on the genome
	p, q := randomTwoInts(genome, rng)
	if p == q {
		return
	}
	// Swap within range
	for i, j := p, q; i < j; i, j = i+1, j-1 {
		genome.Swap(i, j)
	}
	// log.Debugf("After MutInversion: %v", genome)
}

// MutInsertion applies insertion operation on the genome
func MutInsertion(genome eaopt.Slice, rng *rand.Rand) {
	// log.Debugf("Before MutInsertion: %v", genome)
	// Choose two points on the genome
	p, q := randomTwoInts(genome, rng)
	if p == q {
		return
	}

	if rng.Float64() < .5 {
		cq := genome.At(q) // Pop q and insert to p position
		// Move cq to the front and push everyone right
		for i := q; i > p; i-- {
			genome.Set(i, genome.At(i-1))
		}
		genome.Set(p, cq)
	} else {
		cp := genome.At(p)
		// Move cq to the back and push everyone left
		for i := p; i < q; i++ {
			genome.Set(i, genome.At(i+1))
		}
		genome.Set(q, cp)
	}

	// log.Debugf("After MutInsertion: %v", genome)
}

// MutPermute permutes two genes at random n times
func MutPermute(genome eaopt.Slice, rng *rand.Rand) {
	// Nothing to permute
	if genome.Len() <= 1 {
		return
	}
	// Choose two points on the genome
	p, q := randomTwoInts(genome, rng)
	genome.Swap(p, q)
}

// MutSplice splits a genome in 2 and glues the pieces back together in reverse
// order
func MutSplice(genome eaopt.Slice, rng *rand.Rand) {
	var (
		k    = rng.Intn(genome.Len()-1) + 1
		a, b = genome.Split(k)
	)
	genome.Replace(b.Append(a))
}

// Mutate a Tour by applying by inversion or insertion
func (r Tour) Mutate(rng *rand.Rand) {
	rd := rng.Float64()
	if rd < 0.2 {
		MutPermute(r, rng)
	} else if rd < .4 {
		MutSplice(r, rng)
	} else if rd < .7 {
		MutInsertion(r, rng)
	} else {
		MutInversion(r, rng)
	}
}

// Crossover a Tour with another Tour by using Partially Mixed Crossover (PMX).
func (r Tour) Crossover(_ eaopt.Genome, _ *rand.Rand) {
}

// Clone a Tour
func (r Tour) Clone() eaopt.Genome {
	var clone Tour
	clone.Tigs = make([]Tig, r.Len())
	copy(clone.Tigs, r.Tigs)
	clone.M = r.M
	return clone
}

// Shuffle randomly shuffles an integer array using Knuth or Fisher-Yates
func (r Tour) Shuffle(rng *rand.Rand) {
	N := r.Len()
	for i := 0; i < N; i++ {
		// choose index uniformly in [i, N-1]
		j := i + rng.Intn(N-i)
		r.Tigs[j], r.Tigs[i] = r.Tigs[i], r.Tigs[j]
	}
}

// GARun set up the Genetic Algorithm and run it
func (r *CLM) GARun(fwtour *os.File, opt *Optimizer, phase int) (Tour, error) {
	MakeTour := func(rng *rand.Rand) eaopt.Genome {
		c := r.Tour.Clone()
		return c
	}

	ga, err := eaopt.NewDefaultGAConfig().NewGA()
	if err != nil {
		panic(err)
	}

	ga.NPops = 1
	ga.NGenerations = 1000000
	ga.PopSize = uint(opt.NPop)
	ga.Model = eaopt.ModGenerational{
		Selector: eaopt.SelTournament{
			NContestants: 3,
		},
		MutRate: opt.MutProb,
	}
	ga.RNG = opt.rng
	ga.ParallelEval = true

	best := new(float64)
	updated := new(uint)
	*best = -math.MaxFloat64 // Currently best score
	*updated = 0             // Last updated generation

	// Additional bookkeeping per generation
	ga.Callback = func(ga *eaopt.GA) {
		gen := ga.Generations
		currentBest := -ga.HallOfFame[0].Fitness
		if currentBest > *best {
			*best = currentBest
			*updated = gen
		}
		if gen%500 == 0 {
			fmt.Printf("Current iteration GA%d-%d: max_score=%.5f\n",
				phase, gen, currentBest)
			currentBestTour := ga.HallOfFame[0].Genome.(Tour)
			r.printTour(fwtour, currentBestTour, fmt.Sprintf("GA%d-%d-%.5f",
				phase, gen, currentBest))
		}
	}

	// Convergence criteria
	ga.EarlyStop = func(ga *eaopt.GA) bool {
		return ga.Generations-*updated > uint(opt.NGen)
	}

	log.Noticef("GA initialized (npop: %v, ngen: %v, mu: %.2f, rng: %d, break: %d)",
		opt.NPop, opt.NGen, opt.MutProb, opt.Seed, LIMIT)

	err = ga.Minimize(MakeTour)
	r.Tour = ga.HallOfFame[0].Genome.(Tour)
	return r.Tour, err
}
