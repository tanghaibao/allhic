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
	"bytes"
	"strconv"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/xopen"
)

const (
	// LineWidth specifies how many bases to show per line
	LineWidth = 60
	// LargeSequence will notify the writer to send a notification
	LargeSequence = 1000000
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
	lines []AGPLine
}

// NewAGP is the constructor for AGP
func NewAGP(agpfile string) *AGP {
	p := new(AGP)
	return p
}

// Add adds an AGPLine to the collection
func (r *AGP) Add(row string) {
	words := strings.Fields(row)
	var line AGPLine
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
}

// buildFasta builds target FASTA based on info from agpfile
func buildFasta(agpfile string, seqs map[string]*seq.Seq) {
	agp := NewAGP(agpfile)
	file := mustOpen(agpfile)

	log.Noticef("Parse agpfile `%s`", agpfile)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		agp.Add(scanner.Text())
	}

	var buf bytes.Buffer
	outFile := RemoveExt(agpfile) + ".fasta"
	outfh, _ := xopen.Wopen(outFile)
	prevObject := ""
	for _, line := range agp.lines {
		if line.object != prevObject {
			if prevObject != "" {
				writeRecord(prevObject, buf, outfh)
			}
			prevObject = line.object
			buf.Reset()
		}
		if line.isGap {
			buf.Write(bytes.Repeat([]byte("N"), line.gapLength))
		} else {
			if s, ok := seqs[line.componentID]; ok {
				// fmt.Printf("name: %s seq: %s\n", line.componentID, s.SubSeq(1, 10))
				s = s.SubSeq(line.componentBeg, line.componentEnd)
				if line.strand == '-' {
					s.RevComInplace()
				}
				buf.Write(s.Seq)
			} else {
				log.Errorf("Cannot locate %s", line.componentID)
			}
		}
	}
	// Last one
	writeRecord(prevObject, buf, outfh)
	buf.Reset()
	log.Noticef("Assembly FASTA file `%s` built", outFile)
}

// writeRecord writes the FASTA record to the file
func writeRecord(object string, buf bytes.Buffer, outfh *xopen.Writer) {
	record, _ := fastx.NewRecordWithoutValidation(seq.DNA, []byte{}, []byte(object),
		[]byte{}, buf.Bytes())
	size := record.Seq.Length()
	if size > LargeSequence {
		log.Noticef("Write sequence %s (size = %d bp)", record.Name, size)
	}
	record.FormatToWriter(outfh, LineWidth)
	outfh.Flush()
}
