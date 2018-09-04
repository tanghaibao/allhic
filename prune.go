/**
 * Filename: /Users/bao/code/allhic/prune.go
 * Path: /Users/bao/code/allhic
 * Created Date: Wednesday, February 28th 2018, 8:46:25 pm
 * Author: bao
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
)

// Pruner processes the pruning step
type Pruner struct {
	AllelesFile  string
	PairsFile    string
	edges        []ContigPair
	alleleGroups []AlleleGroup
}

// AlleleGroup stores the contig names that are considered allelic
type AlleleGroup []string

// Run calls the pruning steps
// The pruning algorithm is a heuristic method that removes the following pairs:
//
// 1. Allelic, these are directly the pairs of allelic contigs given in the allele table
// 2. Cross-allelic, these are any contigs that connect to the allelic contigs so we only
//    keep the best contig pair
//
// Pruned edges are then annotated as allelic/cross-allelic/ok
func (r *Pruner) Run() {
	r.edges = parseDist(r.PairsFile)
	r.alleleGroups = parseAllelesTable(r.AllelesFile)
	r.pruneAllelic()
	newPairsFile := RemoveExt(r.PairsFile) + ".prune.txt"
	writePairsFile(newPairsFile, r.edges)
}

// pruneAllelic removes the allelic contigs given in the allele table
// we iterate through all the allele groups and mark the pairs that are considered allelic
func (r *Pruner) pruneAllelic() {
	allelicPairs := map[[2]string]bool{}
	for _, alleleGroup := range r.alleleGroups {
		for i := 0; i < len(alleleGroup); i++ {
			for j := i + 1; j < len(alleleGroup); j++ {
				a, b := alleleGroup[i], alleleGroup[j]
				if a > b {
					b, a = a, b
				}
				allelicPairs[[2]string{a, b}] = true
			}
		}
	}

	// Now iterate over all edges and mark
	pruned := 0
	for i, edge := range r.edges {
		pair := [2]string{edge.at, edge.bt}
		if _, ok := allelicPairs[pair]; ok {
			r.edges[i].label = "allelic"
			pruned++
		}
	}
	log.Noticef("Allelic pairs imported: %d, pruned: %s",
		len(allelicPairs), Percentage(pruned, len(allelicPairs)))
}

// parseAllelesTable imports the contig allelic relationship
// File is a tab-separated file that looks like the following:
// Chr10   18902   tig00030660     tig00003333
// Chr10   35071   tig00038687     tig00038686     tig00065419
func parseAllelesTable(allelesFile string) []AlleleGroup {
	log.Noticef("Parse alleles table `%s`", allelesFile)
	fh := mustOpen(allelesFile)
	defer fh.Close()

	reader := bufio.NewReader(fh)

	var data []AlleleGroup
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		words := strings.Split(row, "\t")
		if len(words) <= 3 { // Must have at least 4 fields, i.e. 1 pair
			continue
		}
		alleleGroup := make(AlleleGroup, len(words)-2)
		copy(alleleGroup, words[2:])
		data = append(data, alleleGroup)
	}

	return data
}

// writePairsFile simply writes pruned contig pairs to file
func writePairsFile(pairsFile string, edges []ContigPair) {
	f, _ := os.Create(pairsFile)
	w := bufio.NewWriter(f)
	defer f.Close()
	fmt.Fprintf(w, PairsFileHeader)

	for _, c := range edges {
		fmt.Fprintln(w, c)
	}
	w.Flush()
	log.Noticef("Pruned contig pair analyses written to `%s`", pairsFile)
}
