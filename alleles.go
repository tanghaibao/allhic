/*
 *  alleles.go
 *  allhic
 *
 *  Created by Haibao Tang on 06/28/19
 *  Copyright © 2019 Haibao Tang. All rights reserved.
 */

package allhic

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Alleler is responsible for building the allele table
type Alleler struct {
	PafFile string
}

// Tag represents the additional info in the 12+ columns in the PAF
// file. The type of the tag is dynamically determined
//
// See also:
// https://github.com/lh3/minimap2/blob/master/minimap2.1
//
// The following tags are supported
// Tag	Type	Description
// _
// tp	A	Type of aln: P/primary, S/secondary and I,i/inversion
// cm	i	Number of minimizers on the chain
// s1	i	Chaining score
// s2	i	Chaining score of the best secondary chain
// NM	i	Total number of mismatches and gaps in the alignment
// MD	Z	To generate the ref sequence in the alignment
// AS	i	DP alignment score
// ms	i	DP score of the max scoring segment in the alignment
// nn	i	Number of ambiguous bases in the alignment
// ts	A	Transcript strand (splice mode only)
// cg	Z	CIGAR string (only in PAF)
// cs	Z	Difference string
// dv	f	Approximate per-base sequence divergence
// de	f	Gap-compressed per-base sequence divergence
// rl	i	Length of query regions harboring repetitive seeds
type Tag = interface{}

// PAFRecord holds one line in the PAF file
// The file spec:
// https://github.com/lh3/miniasm/blob/master/PAF.md
type PAFRecord struct {
	query           string         // Query sequence name
	queryLength     int            // Query sequence length
	queryStart      int            // Query start (0-based)
	queryEnd        int            // Query end (0-based)
	relativeStrand  byte           // `+' if query and target on the same strand; `-' if opposite
	target          string         // Target sequence name
	targetLength    int            // Target sequence length
	targetStart     int            // Target start on original strand (0-based)
	targetEnd       int            // Target end on original strand (0-based)
	nMatches        int            // Number of matching bases in the mapping
	alignmentLength int            // Number bases, including gaps, in the mapping
	mappingQuality  uint8          // Mapping quality (0-255 with 255 for missing)
	tags            map[string]Tag // Tags, e.g. tp, cm etc.
}

// PAF parses the PAF file into a set of records
type PAF struct {
	PafFile string      // File path of the paf
	records []PAFRecord // List of PAF records
}

// ParseRecords collects all records in memory
func (r *PAF) ParseRecords() {
	r.records = []PAFRecord{}
	fh := mustOpen(r.PafFile)

	log.Noticef("Parse paffile `%s`", r.PafFile)
	reader := bufio.NewReader(fh)
	var rec PAFRecord

	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		words := strings.Split(row, "\t")
		if len(words) < 12 || len(words[4]) < 1 {
			continue
		}

		// Parse the first 12 columns
		rec.query = words[0]
		rec.queryLength, _ = strconv.Atoi(words[1])
		rec.queryStart, _ = strconv.Atoi(words[2])
		rec.queryEnd, _ = strconv.Atoi(words[3])
		rec.relativeStrand = words[4][0]
		rec.target = words[5]
		rec.targetLength, _ = strconv.Atoi(words[6])
		rec.targetStart, _ = strconv.Atoi(words[7])
		rec.targetEnd, _ = strconv.Atoi(words[8])
		rec.nMatches, _ = strconv.Atoi(words[9])
		rec.alignmentLength, _ = strconv.Atoi(words[10])
		mappingQuality, _ := strconv.Atoi(words[11])
		rec.mappingQuality = uint8(mappingQuality)
		rec.tags = map[string]Tag{}
		var tag Tag

		// Parse columns 12+
		for i := 12; i < len(words); i++ {
			tokens := strings.Split(words[i], ":")
			if len(tokens) < 3 {
				continue
			}
			tagName := tokens[0]
			value := tokens[2]
			switch tokens[1] {
			case "i":
				tag, _ = strconv.Atoi(value)
			case "f":
				tag, _ = strconv.ParseFloat(value, 32)
			default:
				tag = value
			}
			rec.tags[tagName] = tag
		}

		r.records = append(r.records, rec)
	}
}

// Run kicks off the Alleler
func (r *Alleler) Run() {
	paf := PAF{PafFile: r.PafFile}
	paf.ParseRecords()
	fmt.Println(paf.records)
	log.Notice("Success")
}
