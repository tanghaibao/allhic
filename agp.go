/**
 * Filename: /Users/bao/code/allhic/agp.go
 * Path: /Users/bao/code/allhic
 * Created Date: Monday, February 26th 2018, 8:30:12 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

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
	componentID  string
	componentBeg int
	componentEnd int
}

// AGP is a collection of AGPLines
type AGP struct {
	lines []*AGPLine
}

// NewAGP is the constructor for AGP
func NewAGP(agpfile string) *AGP {
	p := new(AGP)
	return p
}

// Add adds an AGPLine to the collection
func (r *AGP) Add(row string) {
	words := strings.Fields(row)
	line := new(AGPLine)
	line.object = words[0]
	line.objectBeg, _ = strconv.Atoi(words[1])
	line.objectEnd, _ = strconv.Atoi(words[2])
	line.partNumber, _ = strconv.Atoi(words[3])
	line.componentType = words[4][0]
	line.isGap = (line.componentType == 'N' || line.componentType == 'U')
	if line.isGap {
		line.gapLength, _ = strconv.Atoi(words[5])
		line.gapType = words[6]
		line.linkage = words[7]
		line.linkageEvidence = words[8]
	} else {
		line.componentID = words[5]
		line.componentBeg, _ = strconv.Atoi(words[6])
		line.componentEnd, _ = strconv.Atoi(words[7])
		line.strand = words[8][0]
	}
	r.lines = append(r.lines, line)
	fmt.Println(*line)
}

// BuildFasta builds target FASTA based on info from agpfile
func BuildFasta(agpfile, fastafile string) {
	agp := NewAGP(agpfile)
	file, _ := os.Open(agpfile)

	log.Noticef("Parse agpfile `%s`", agpfile)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		agp.Add(scanner.Text())
	}
}
