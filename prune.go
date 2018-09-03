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
	"strings"
)

// Pruner processes the pruning step
type Pruner struct {
	AllelesFile string
	PairsFile   string
}

// AlleleGroup stores the contig names that are considered allelic
type AlleleGroup []string

// Run calls the pruning steps
func (r *Pruner) Run() {
	edges := parseDist(r.PairsFile)
	alleleGroups := parseAllelesTable(r.AllelesFile)
	fmt.Println(edges[0])
	fmt.Println(alleleGroups[0])
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
