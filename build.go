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
	componentType byte
	isGap         bool
	strand        byte
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
	strand        byte
}

// OO describes a scaffolding experiment and contains an array of OOLine
type OO struct {
	sizes   map[string]int
	entries []OOLine
}

// IsNewerFile checks if file a is newer than file b
func IsNewerFile(a, b string) bool {
	af, aerr := os.Stat(a)
	bf, berr := os.Stat(b)
	if os.IsNotExist(aerr) || os.IsNotExist(berr) {
		return false
	}
	am := af.ModTime()
	bm := bf.ModTime()
	return am.Sub(bm) > 0
}

// GetFastaSizes returns a dictionary of contig sizes
func (r *OO) GetFastaSizes(fastafile string) {
	log.Noticef("Parse FASTA file `%s`", fastafile)
	r.sizes = make(map[string]int)
	faifile := fastafile + ".fai"
	// Check if the .fai file is outdated
	if !IsNewerFile(faifile, fastafile) {
		os.Remove(faifile)
	}

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
func (r *OO) Add(scaffold, ctg string, ctgsize int, strand byte) {
	o := OOLine{scaffold, ctg, ctgsize, strand}
	r.entries = append(r.entries, o)
}

// WriteAGP converts the simplistic OOLine into AGP format
func (r *Builder) WriteAGP(oo *OO, filename string) {
	gapSize := 100
	gapType := "scaffold"
	linkage := "yes"
	evidence := "map"
	log.Noticef("Gapsize = %d, Gaptype = %s, Linkage = %s", gapSize, gapType, linkage)

	prevObject := ""
	objectBeg := 1
	objectEnd := 1
	partNumber := 0
	componentType := 'W'
	f, _ := os.Create(filename)
	defer f.Close()
	w := bufio.NewWriter(f)
	components := 0

	// Write AGP for each object group
	for _, line := range oo.entries {
		if line.id != prevObject {
			prevObject = line.id
			objectBeg = 1
			partNumber = 0
		}
		if partNumber > 0 && gapSize > 0 {
			if gapSize == 100 {
				componentType = 'U'
			} else {
				componentType = 'N'
			}
			objectEnd = objectBeg + gapSize - 1
			partNumber++
			fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%c\t%d\t%s\t%s\t%s\n",
				line.id, objectBeg, objectEnd, partNumber,
				componentType, gapSize, gapType, linkage, evidence)
			objectBeg += gapSize
		}
		objectEnd = objectBeg + line.componentSize - 1
		partNumber++
		fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%c\n",
			line.id, objectBeg, objectEnd, partNumber,
			'W', line.componentID, 1, line.componentSize, line.strand)
		objectBeg += line.componentSize
		components++
	}
	w.Flush()
	log.Noticef("A total of %d tigs written to `%s`", components, filename)
}

// Run kicks off the Builder
func (r *Builder) Run() {
	oo := r.ReadFiles()
	agpfile := RemoveExt(r.Tourfile) + ".agp"
	r.WriteAGP(oo, agpfile)
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
	var (
		name   string
		strand byte
	)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		if words[0][0] == '>' {
			name = words[0][1:]
			continue
		}
		for _, tig := range words {
			at, ao := tig[:len(tig)-1], tig[len(tig)-1]
			if ao == '+' || ao == '-' || ao == '?' {
				tig, strand = at, ao
			} else {
				strand = '?'
			}
			r.Add(name, tig, r.sizes[tig], strand)
		}
	}
	fmt.Println(r)
}
