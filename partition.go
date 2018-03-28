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
	r.CountLinks()
	// fmt.Println(counts)
	log.Notice("Success")
}

// CountLinks provides the method to count the links
func (r *Partitioner) CountLinks() [][]int {
	fh, _ := os.Open(r.Bamfile)
	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	tigToIdx := make(map[string]int)
	refs := br.Header().Refs()
	for i, ref := range refs {
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
		C[ai][bi]++
		C[bi][ai]++
	}

	// Write the edge list
	for i := 0; i < N; i++ {
		for j := i + 1; j < N; j++ {
			if C[i][j] != 0 {
				fmt.Printf("%s\t%s\t%d\n", refs[i].Name(), refs[j].Name(), C[i][j])
			}
		}
	}
	return C
}
