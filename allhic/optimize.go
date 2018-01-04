/**
 * Filename: /Users/bao/code/allhic/allhic/optimize.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Tuesday, January 2nd 2018, 10:00:33 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import "fmt"

// Optimizer runs the order-and-orientation procedure, given a clmfile
type Optimizer struct {
	Clmfile string
}

// Run kicks off the Optimizer
func (r *Optimizer) Run() {
	clm := InitCLMFile(r.Clmfile)
	tour := clm.Activate()
	fmt.Println(tour)

	//M := clm.M()
	//fmt.Println(M)
}
