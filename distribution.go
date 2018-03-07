/**
 * Filename: /Users/bao/code/allhic/distribution.go
 * Path: /Users/bao/code/allhic
 * Created Date: Wednesday, March 7th 2018, 1:56:45 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import "fmt"

// Distribution processes the distribution step
type Distribution struct {
	Bamfile string
}

// Run calls the distribution steps
func (r *Distribution) Run() {
	fmt.Println(r.Bamfile)
}
