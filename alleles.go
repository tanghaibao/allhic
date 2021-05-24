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
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// Alleler is responsible for building the allele table
type Alleler struct {
	PafFile  string       // ex. "genome.paf"
	ReFile   string       // ex. "genome.counts_GATC.txt"
	Paf      PAFFile      // The PAF data
	ReCounts RECountsFile // The RE data
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
	Query           string         // Query sequence name
	QueryLength     int            // Query sequence length
	QueryStart      int            // Query start (0-based)
	QueryEnd        int            // Query end (0-based)
	RelativeStrand  byte           // `+' if query and target on the same strand; `-' if opposite
	Target          string         // Target sequence name
	TargetLength    int            // Target sequence length
	TargetStart     int            // Target start on original strand (0-based)
	TargetEnd       int            // Target end on original strand (0-based)
	NumMatches      int            // Number of matching bases in the mapping
	AlignmentLength int            // Number bases, including gaps, in the mapping
	MappingQuality  uint8          // Mapping quality (0-255 with 255 for missing)
	Tags            map[string]Tag // Tags, e.g. tp, cm etc.
}

// PAFFile parses the PAF file into a set of records
type PAFFile struct {
	PafFile string      // File path of the paf
	Records []PAFRecord // List of PAF records
}

// ParseRecords collects all records in memory
func (r *PAFFile) ParseRecords() error {
	r.Records = []PAFRecord{}
	fh, err := os.Open(r.PafFile)
	if err != nil {
		return err
	}

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
		rec.Query = words[0]
		rec.QueryLength, _ = strconv.Atoi(words[1])
		rec.QueryStart, _ = strconv.Atoi(words[2])
		rec.QueryEnd, _ = strconv.Atoi(words[3])
		rec.RelativeStrand = words[4][0]
		rec.Target = words[5]
		rec.TargetLength, _ = strconv.Atoi(words[6])
		rec.TargetStart, _ = strconv.Atoi(words[7])
		rec.TargetEnd, _ = strconv.Atoi(words[8])
		rec.NumMatches, _ = strconv.Atoi(words[9])
		rec.AlignmentLength, _ = strconv.Atoi(words[10])
		mappingQuality, _ := strconv.Atoi(words[11])
		rec.MappingQuality = uint8(mappingQuality)
		rec.Tags = map[string]Tag{}
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
			rec.Tags[tagName] = tag
		}

		r.Records = append(r.Records, rec)
	}
	return nil
}

// extractAllelicPairs collects Extract allelic pairs
func (r *Alleler) extractAllelicPairs() {
	// Sort the contigs by sizes, starting from shortest
	sort.Slice(r.ReCounts.Records, func(i, j int) bool {
		return r.ReCounts.Records[i].Length < r.ReCounts.Records[j].Length
	})
	// Find significant matches of small-big allelic contig pairs
}

// Run kicks off the Alleler
func (r *Alleler) Run() error {
	r.Paf = PAFFile{PafFile: r.PafFile}
	err := r.Paf.ParseRecords()
	if err != nil {
		return err
	}
	r.ReCounts = RECountsFile{Filename: r.ReFile}
	err = r.ReCounts.ParseRecords()
	if err != nil {
		return err
	}
	r.extractAllelicPairs()
	log.Notice("Success")
	return nil
}
