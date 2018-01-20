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
	"io"
	"os"

	"github.com/biogo/hts/bam"
)

// Partitioner converts the bamfile into a matrix of link counts
type Partitioner struct {
	Bamfile string
}

// Run is the main function body of partition
func (r *Partitioner) Run() {
	counts := r.CountLinks()
	fmt.Println(counts)
	log.Notice("Success")
}

// CountLinks provides the method to count the links
func (r *Partitioner) CountLinks() [][]int {
	fh, _ := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	tigToIdx := make(map[string]int)
	for i, ref := range br.Header().Refs() {
		tigToIdx[ref.Name()] = i
	}

	N := len(tigToIdx)
	C := Make2DSlice(N, N)
	log.Noticef("Initiating matrix of size %d x %d", N, N)
	// Collect all links in a 2D matrix
	for {
		rec, err := br.Read()
		if err != nil {
			if err != io.EOF {
				log.Error(err)
			}
			break
		}
		ai := tigToIdx[rec.Ref.Name()]
		bi := tigToIdx[rec.MateRef.Name()]
		// fmt.Println(rec.Ref.Name(), rec.MateRef.Name())
		C[ai][bi]++
		C[bi][ai]++
	}
	return C
}
