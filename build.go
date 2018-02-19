/**
 * Filename: /Users/htang/code/allhic/build.go
 * Path: /Users/htang/code/allhic
 * Created Date: Saturday, January 27th 2018, 10:21:08 pm
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"fmt"

	"github.com/shenwei356/bio/seqio/fai"
)

// Builder reconstructs the genome release AGP and FASTA files
type Builder struct {
	Tourfile  string
	Fastafile string
}

// AGPLine is a line in the AGP file
type AGPLine struct {
	object        string
	objectBeg     int
	objectEnd     int
	partNumber    int
	componentType rune
	isGap         bool
	orientation   string
	// As a gap
	gapLength       int
	gapType         string
	linkage         string
	linkageEvidence string
	// As a sequence chunk
	componentID  int
	componentBeg int
	componentEnd int
}

// OOLine describes a simple contig entry in a scaffolding experiment
type OOLine struct {
	id            string
	componentID   string
	componentSize int
	strand        string
}

// OO describes a scaffolding experiment and contains an array of OOLine
type OO struct {
	contigs []string
	sizes   map[string]int
	entries []OOLine
}

// GetFastaSizes returns a dictionary of contig sizes
func GetFastaSizes(fastafile string) map[string]int {
	sizes := make(map[string]int)
	faidx, _ := fai.New(fastafile)
	defer faidx.Close()

	for name, rec := range faidx.Index {
		sizes[name] = rec.Length
	}

	return sizes
}

// NewOO initializes OO object
func NewOO(fastafile string) *OO {
	p := new(OO)
	p.sizes = GetFastaSizes(fastafile)
	return p
}

// Add instantiates a new OOLine object and add to the array in OO
func (r *OO) Add(scaffold, ctg string, ctgsize int, strand string) {
	o := OOLine{scaffold, ctg, ctgsize, strand}
	r.entries = append(r.entries, o)
}

// WriteAGP converts the simplistic OOLine into AGP format
func (r *OO) WriteAGP() {
	gapSize := 100
	gapType := "scaffold"
	log.Noticef("Gapsize = %d, Gaptype = %s", gapSize, gapType)
	fmt.Println(r.sizes)
}

// Run kicks off the Builder
func (r *Builder) Run() {
	fmt.Println(r.Tourfile, r.Fastafile)
	oo := NewOO(r.Fastafile)
	oo.WriteAGP()
}
