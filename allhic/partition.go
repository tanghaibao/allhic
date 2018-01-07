/**
 * Filename: /Users/htang/code/allhic/allhic/partition.go
 * Path: /Users/htang/code/allhic/allhic
 * Created Date: Wednesday, January 3rd 2018, 11:21:45 am
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"fmt"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// Partitioner converts the bamfile into a matrix of link counts
type Partitioner struct {
	Bamfile string
}

// Run is the main function body of partition
func (r *Partitioner) Run() {
	r.CountLinks()
	log.Notice("Success")
}

// CountLinks provides the method to count the links
func (r *Partitioner) CountLinks() {
	fh, _ := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	for _, rg := range br.Header().RGs() {
		fmt.Println(rg.Get(sam.Tag([2]byte{'S', 'M'})))
	}
}
