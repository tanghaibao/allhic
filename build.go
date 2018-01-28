/**
 * Filename: /Users/htang/code/allhic/build.go
 * Path: /Users/htang/code/allhic
 * Created Date: Saturday, January 27th 2018, 10:21:08 pm
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import "fmt"

// Builder reconstructs the genome release AGP and FASTA files
type Builder struct {
	Tourfile  string
	Fastafile string
}

// Run kicks off the Builder
func (r *Builder) Run() {
	fmt.Println(r.Tourfile, r.Fastafile)
}
