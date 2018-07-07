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
	"io"
	"os"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
)

// Builder reconstructs the genome release AGP and FASTA files
type Builder struct {
	Tourfile  string
	Fastafile string
	AGPfile   string
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
	seqs    map[string]*seq.Seq
	entries []OOLine
}

// getFastaSizes returns a dictionary of contig sizes
func (r *OO) getFastaSizes(fastafile string) {
	log.Noticef("Parse FASTA file `%s`", fastafile)

	reader, _ := fastx.NewDefaultReader(fastafile)
	seq.ValidateSeq = false
	r.seqs = map[string]*seq.Seq{}
	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}

		name := string(rec.Name)
		r.seqs[name] = rec.Seq.Clone()
	}
}

// readFiles initializes OO object
func (r *Builder) readFiles() *OO {
	oo := new(OO)
	oo.getFastaSizes(r.Fastafile)
	oo.parseLastTour(r.Tourfile)

	return oo
}

// Add instantiates a new OOLine object and add to the array in OO
func (r *OO) Add(scaffold, ctg string, ctgsize int, strand byte) {
	o := OOLine{scaffold, ctg, ctgsize, strand}
	r.entries = append(r.entries, o)
}

// writeAGP converts the simplistic OOLine into AGP format
func (r *Builder) writeAGP(oo *OO) {
	filename := r.AGPfile
	gapSize := 100
	gapType := "scaffold"
	linkage := "yes"
	evidence := "map"
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
	r.AGPfile = RemoveExt(r.Tourfile) + ".agp"
	r.Build(r.readFiles())
	log.Notice("Success")
}

// Build constructs molecule using component FASTA sequence
func (r *Builder) Build(oo *OO) {
	r.writeAGP(oo)
	buildFasta(r.AGPfile, oo.seqs)
}

// parseLastTour reads tour from file
//
// A tour file has the following format:
// > name
// contig1+ contig2- contig3?
func (r *OO) parseLastTour(tourfile string) {
	words := parseTourFile(tourfile)
	var strand byte
	for _, tig := range words {
		at, ao := tig[:len(tig)-1], tig[len(tig)-1]
		if ao == '+' || ao == '-' || ao == '?' {
			tig, strand = at, ao
		} else {
			strand = '?'
		}
		r.Add("FINAL", tig, r.seqs[tig].Length(), strand)
	}
}

// ParseAllTours reads tour from file
//
// A tour file has the following format:
// > name
// contig1+ contig2- contig3?
func (r *OO) ParseAllTours(tourfile string) {
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
			r.Add(name, tig, r.seqs[tig].Length(), strand)
		}
	}
}
