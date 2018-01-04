/**
 * Filename: /Users/bao/code/allhic/allhic/evaluate.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Wednesday, January 3rd 2018, 9:40:36 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import "fmt"

const limit = 10000000

// EvaluateM calculates a score for the current tour
func EvaluateM(tour, tourSizes []int, M [][]int) float64 {
	size := len(tour)
	mid := make([]float64, size)
	cumSum := 0.0
	for i := range tour {
		tsize := float64(tourSizes[i])
		mid[i] = cumSum + tsize/2
		cumSum += tsize
	}
	fmt.Println(mid)

	// Now add up all the pairwise scores
	s := 0.0
	for i := 0; i < size; i++ {
		for j := i + 1; j < size; j++ {
			nlinks := M[tour[i]][tour[j]]
			dist := mid[j] - mid[i]
			if dist > limit {
				break
			}
			s += float64(nlinks) / dist
		}
	}
	return s
}
