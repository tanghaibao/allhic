/**
 * Filename: /Users/bao/code/allhic/allhic/optimize.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Tuesday, January 2nd 2018, 10:00:33 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"fmt"
	"math/rand"
)

// Optimizer runs the order-and-orientation procedure, given a clmfile
type Optimizer struct {
	Clmfile string
}

// Shuffle randomly shuffles an integer array using Knuth or Fisher-Yates
func Shuffle(tour Tour) {
	N := len(tour.Tigs)
	for i := 0; i < N; i++ {
		// choose index uniformly in [i, N-1]
		r := i + rand.Intn(N-i)
		tour.Tigs[r], tour.Tigs[i] = tour.Tigs[i], tour.Tigs[r]
	}
}

// Run kicks off the Optimizer
func (r *Optimizer) Run() {
	clm := InitCLMFile(r.Clmfile)
	tour := clm.Activate()

	// fmt.Println(tour)
	tourScore := EvaluateM(tour)
	fmt.Println(tourScore)

	// var test []string
	// for i := range tour {
	// 	test = append(test, clm.activeTigs[i])
	// }
	// fmt.Println(test)

	Shuffle(tour)
	// fmt.Println(tour)
	tourScore = EvaluateM(tour)
	fmt.Println(tourScore)
}
