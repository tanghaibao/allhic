/**
 * Filename: /Users/bao/code/allhic/prune.go
 * Path: /Users/bao/code/allhic
 * Created Date: Wednesday, February 28th 2018, 8:46:25 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import "fmt"

// Pruner processes the pruning step
type Pruner struct {
	Bamfile string
}

// Run calls the pruning steps
func (r *Pruner) Run() {
	fmt.Println(r.Bamfile)
}
