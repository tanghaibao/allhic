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
	"bufio"
	"encoding/csv"
	"io"
	"os"
	"strconv"
)

// Partitioner converts the bamfile into a matrix of link counts
type Partitioner struct {
	Distfile string
}

// Run is the main function body of partition
func (r *Partitioner) Run() {
	r.ParseDist()
	log.Notice("Success")
}

// ParseDist imports the edges of the contig linkage graph
func (r *Partitioner) ParseDist() {
	edges := ParseDistLines(r.Distfile)
	goodEdges := FilterEdges(edges)
	// fmt.Println(goodEdges)
	log.Noticef("Edge filtering keeps %s edges",
		Percentage(len(goodEdges), len(edges)))
}

// FilterEdges implements rules to keep edges between close contigs and remove distant or weak contig pairs
func FilterEdges(edges []ContigPair) []ContigPair {
	var goodEdges []ContigPair

	for _, e := range edges {
		if e.mleDistance >= EffLinkDist {
			continue
		}
		goodEdges = append(goodEdges, e)
	}

	return goodEdges
}

// ParseDistLines imports the edges of the contig into a slice of DistLine
// DistLine stores the data structure of the Dist file
// #Contig1        Contig2 Length1 Length2 LDE1    LDE2    LDE     ObservedLinks   ExpectedLinksIfAdjacent MLEdistance
// jpcChr1.ctg199  jpcChr1.ctg257  124567  274565  0.3195  2.0838  1.1607  2       27.4    1617125
// idcChr1.ctg353  idcChr1.ctg382  143105  270892  2.1577  1.0544  1.3505  2       34.2    2190000
func ParseDistLines(distfile string) []ContigPair {
	var edges []ContigPair

	fh, _ := os.Open(distfile)
	defer fh.Close()
	r := csv.NewReader(bufio.NewReader(fh))
	r.Comma = '\t'
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		at, bt := rec[0], rec[1]
		L1, _ := strconv.Atoi(rec[2])
		L2, _ := strconv.Atoi(rec[3])
		lde1, _ := strconv.ParseFloat(rec[4], 64)
		lde2, _ := strconv.ParseFloat(rec[5], 64)
		localLDE, _ := strconv.ParseFloat(rec[6], 64)
		nObservedLinks, _ := strconv.Atoi(rec[7])
		nExpectedLinks, _ := strconv.ParseFloat(rec[8], 64)
		mleDistance, _ := strconv.Atoi(rec[9])

		cp := ContigPair{
			at: at, bt: bt,
			L1: L1, L2: L2,
			lde1: lde1, lde2: lde2, localLDE: localLDE,
			nObservedLinks: nObservedLinks, nExpectedLinks: nExpectedLinks,
			mleDistance: mleDistance}

		edges = append(edges, cp)
	}

	return edges
}
