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
	"bufio"
	"fmt"
	"os"
	"strings"

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
	name         string
	contigs      []string
	orientations []byte
	sizes        map[string]int
	entries      []OOLine
}

// GetFastaSizes returns a dictionary of contig sizes
func (r *OO) GetFastaSizes(fastafile string) {
	log.Noticef("Parse FASTA file `%s`", fastafile)
	r.sizes = make(map[string]int)
	faidx, _ := fai.New(fastafile)
	defer faidx.Close()

	for name, rec := range faidx.Index {
		r.sizes[name] = rec.Length
	}
}

// ReadFiles initializes OO object
func (r *Builder) ReadFiles() *OO {
	oo := new(OO)
	oo.GetFastaSizes(r.Fastafile)
	oo.ParseTour(r.Tourfile)

	return oo
}

// Add instantiates a new OOLine object and add to the array in OO
func (r *OO) Add(scaffold, ctg string, ctgsize int, strand string) {
	o := OOLine{scaffold, ctg, ctgsize, strand}
	r.entries = append(r.entries, o)
}

// WriteAGP converts the simplistic OOLine into AGP format
func (r *Builder) WriteAGP(oo *OO) {
	gapSize := 100
	gapType := "scaffold"
	log.Noticef("Gapsize = %d, Gaptype = %s", gapSize, gapType)
	fmt.Println(oo.sizes)
}

// Run kicks off the Builder
func (r *Builder) Run() {
	oo := r.ReadFiles()
	r.WriteAGP(oo)
}

// ParseTour reads tour from file
//
// A tour file has the following format:
// > name
// contig1+ contig2- contig3?
func (r *OO) ParseTour(tourfile string) {
	log.Noticef("Parse tourfile `%s`", tourfile)

	file, _ := os.Open(tourfile)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		if words[0][0] == '>' {
			r.name = words[0][1:]
		} else {
			for _, tig := range words {
				at, ao := tig[:len(tig)-1], tig[len(tig)-1]
				if ao == '+' || ao == '-' || ao == '?' {
					r.contigs = append(r.contigs, at)
					r.orientations = append(r.orientations, ao)
				} else {
					r.contigs = append(r.contigs, tig)
					r.orientations = append(r.orientations, '?')
				}
			}
		}
	}
	fmt.Println(r)
}
