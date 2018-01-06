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

// CountLinks provides the method to count the links
func (r Partitioner) CountLinks() {
	fh, err := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	if err != nil {
		panic(err)
	}

	br, err := bam.NewReader(fh, 0)
	if err != nil {
		panic(err)
	}
	defer br.Close()

	for _, rg := range br.Header().RGs() {
		fmt.Println(rg.Get(sam.Tag([2]byte{'S', 'M'})))
	}
}
