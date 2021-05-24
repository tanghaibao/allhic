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

	hungarianAlgorithm "github.com/oddg/hungarian-algorithm"
)

// Pruner processes the pruning step
type Pruner struct {
	AllelesFile  string
	PairsFile    string
	edges        []ContigPair
	alleleGroups []AlleleGroup
}

// ContigAB is used to get a pair of contigs
type ContigAB [2]string

// AlleleGroup stores the contig names that are considered allelic
type AlleleGroup []string

// CtgAlleleGroupPair stores a pair of the contig and the alleleGroup it resides in
type CtgAlleleGroupPair struct {
	ctg     string
	groupID int
}

// Run calls the pruning steps
// The pruning algorithm is a heuristic method that removes the following pairs:
//
// 1. Allelic, these are directly the pairs of allelic contigs given in the allele table
// 2. Cross-allelic, these are any contigs that connect to the allelic contigs so we only
//    keep the best contig pair
//
// Pruned edges are then annotated as allelic/cross-allelic/ok
func (r *Pruner) Run() error {
	edges, err := parseDist(r.PairsFile)
	if err != nil {
		return err
	}
	r.edges = edges
	alleleGroups, err := parseAllelesFile(r.AllelesFile)
	if err != nil {
		return err
	}
	r.alleleGroups = alleleGroups
	r.pruneAllelic()
	r.pruneCrossAllelicBipartiteMatching()
	// r.pruneCrossAllelic()
	newPairsFile := RemoveExt(r.PairsFile) + ".prune.txt"
	return writePairsFile(newPairsFile, r.edges)
}

// pruneAllelic removes the allelic contigs given in the allele table
// we iterate through all the allele groups and mark the pairs that are considered allelic
func (r *Pruner) pruneAllelic() {
	// Find all blacklisted allelic pairs
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
	pruned, prunedLinks := 0, 0
	total, totalLinks := 0, 0
	for i, edge := range r.edges {
		pair := [2]string{edge.at, edge.bt}
		if _, ok := allelicPairs[pair]; ok {
			r.edges[i].label = "allelic"
			pruned++
			prunedLinks += edge.nObservedLinks
		}
		total++
		totalLinks += edge.nObservedLinks
	}
	log.Noticef("Allelic pairs imported: %d, pruned: %s, prunedLinks: %s",
		len(allelicPairs), Percentage(pruned, total), Percentage(prunedLinks, totalLinks))
}

// pruneCrossAllelicBipartiteMatching is a heuristic that tests whether an edge
// is weak based on maximum weight bipartite matching. For example:
//
// a ===== 100 ===== b
//
// (a-d) 20  (b-c) 120
//
// c ===== 100 ===== d
// we still favor the partition where ab + cd = 200 > ad + bc = 140
func (r *Pruner) pruneCrossAllelicBipartiteMatching() {
	ctgToAlleleGroup := r.getCtgToAlleleGroup()

	// Store scores for contig pairs
	ctgPairScores := map[ContigAB]int{}
	for _, edge := range r.edges { // First pass collects all scores
		if edge.label != "ok" { // We skip the allelic pairs since these are already removed
			continue
		}
		ctgPairScores[ContigAB{edge.at, edge.bt}] = edge.nObservedLinks
		ctgPairScores[ContigAB{edge.bt, edge.at}] = edge.nObservedLinks
	}

	// Now iterate over all edges and mark
	pruned, prunedLinks := 0, 0
	total, totalLinks := 0, 0
	for i, edge := range r.edges {
		if edge.label != "ok" {
			continue
		}
		if !r.isStrongEdgeInBipartiteMatchingGroups(&edge, ctgToAlleleGroup, ctgPairScores) {
			r.edges[i].label = edge.label
			pruned++
			prunedLinks += edge.nObservedLinks
		}
		total++
		totalLinks += edge.nObservedLinks
	}
	log.Noticef("Cross-allelic pairs pruned: %s, prunedLinks: %s",
		Percentage(pruned, total), Percentage(prunedLinks, totalLinks))
}

// isStrongEdgeInBipartiteMatching determines if the edge being considered is
// used in bipartite matching between two allele groups on either side of this
// edge, since a-b can both be within a number of AlleleGroups. We need to check
// each pair one by one.
func (r *Pruner) isStrongEdgeInBipartiteMatchingGroups(edge *ContigPair, ctgToAlleleGroup map[string][]int, ctgPairScores map[ContigAB]int) bool {
	ag, aok := ctgToAlleleGroup[edge.at]
	bg, bok := ctgToAlleleGroup[edge.bt]
	if !aok || !bok {
		return true
	}
	for _, ai := range ag {
		for _, bi := range bg {
			aGroup := r.alleleGroups[ai]
			bGroup := r.alleleGroups[bi]
			if !r.isStrongEdgeInBipartiteMatching(edge, aGroup, bGroup, ctgPairScores) {
				edge.label = fmt.Sprintf("cross-allelic(%s|%s)", strings.Join(aGroup, ","), strings.Join(bGroup, ","))
				return false
			}
		}
	}
	return true
}

// isStrongEdgeInBipartiteMatching determines if the edge being considered is
// used in bipartite matching between two allele groups on either side of this
// edge. Note that this function is called by
// isStrongEdgeInBipartiteMatchingGroups(), and only operates on a single pair
// of AlleleGroups.
func (r *Pruner) isStrongEdgeInBipartiteMatching(edge *ContigPair, aGroup AlleleGroup, bGroup AlleleGroup, ctgPairScores map[ContigAB]int) bool {
	// Build a square matrix that contain matching scores
	aN := len(aGroup)
	bN := len(bGroup)
	N := max(aN, bN)
	S := Make2DSlice(N, N)
	ti, tj := -1, -1
	// Populate the entries
	for i, at := range aGroup {
		if at == edge.at {
			ti = i
		}
		for j, bt := range bGroup {
			if bt == edge.bt {
				tj = j
			}
			ctgPair := ContigAB{at, bt}
			if score, ok := ctgPairScores[ctgPair]; ok {
				S[i][j] = score
			}
		}
	}
	// Solve the matching problem using Hungarian algorithm
	solution := maxBipartiteMatchingWithWeights(S)
	ans := solution[ti] == tj
	// fmt.Println(edge.at, edge.bt, aGroup, bGroup, S, solution, ans)
	return ans
}

// maxBipartiteMatchingWithWeights calculates the bipartite matching using the
// weights, wraps hungarianAlgorithm() which minimizes the costs, so we need to
// transform from weights to costs
func maxBipartiteMatchingWithWeights(weights [][]int) []int {
	maxCell := 0
	// Get the max value of the matrix
	for _, row := range weights {
		for _, cell := range row {
			maxCell = max(maxCell, cell)
		}
	}
	N := len(weights)
	costs := Make2DSlice(N, N)
	// Subtract the weights from the max to get costs
	for i, row := range weights {
		for j, cell := range row {
			costs[i][j] = maxCell - cell
		}
	}
	// By default, hungarianAlgorithm works on costs
	solution, _ := hungarianAlgorithm.Solve(costs)
	return solution
}

// getCtgToAlleleGroup returns contig to List of alleleGroups
func (r *Pruner) getCtgToAlleleGroup() map[string][]int {
	// Store contig to list of alleleGroups since each contig can be in different alleleGroups
	ctgToAlleleGroup := map[string][]int{}
	for groupID, alleleGroup := range r.alleleGroups {
		for _, ctg := range alleleGroup {
			if gg, ok := ctgToAlleleGroup[ctg]; ok {
				gg = append(gg, groupID)
			} else {
				ctgToAlleleGroup[ctg] = []int{groupID}
			}
		}
	}
	return ctgToAlleleGroup
}

// pruneCrossAllelic removes contigs that link to multiple allelic contigs and we choose
// to keep a single best link, i.e. for each allele group, we find the best link to each
// contig and retain. Note that since contigs may be in different allele groups, the order
// of removal may affect end results.
func (r *Pruner) pruneCrossAllelic() {
	ctgToAlleleGroup := r.getCtgToAlleleGroup()

	// Store the best match of each contig to an allele group
	scores := map[CtgAlleleGroupPair]int{} // (ctg, alleleGroupID) => score
	for _, edge := range r.edges {
		if edge.label != "ok" { // We skip the allelic pairs since these are already removed
			continue
		}
		updateScore(edge.at, edge.bt, edge.nObservedLinks, ctgToAlleleGroup, scores)
		updateScore(edge.bt, edge.at, edge.nObservedLinks, ctgToAlleleGroup, scores)
	}

	// Now iterate over all edges and mark
	pruned, prunedLinks := 0, 0
	total, totalLinks := 0, 0
	for i, edge := range r.edges {
		aBestScore := getScore(edge.at, edge.bt, ctgToAlleleGroup, scores)
		bBestScore := getScore(edge.bt, edge.at, ctgToAlleleGroup, scores)
		if edge.nObservedLinks < aBestScore && edge.nObservedLinks < bBestScore {
			r.edges[i].label = fmt.Sprintf("cross-allelic(%d|%d)", aBestScore, bBestScore)
			pruned++
			prunedLinks += edge.nObservedLinks
		}
		total++
		totalLinks += edge.nObservedLinks
	}
	log.Noticef("Cross-allelic pairs pruned: %s, prunedLinks: %s",
		Percentage(pruned, total), Percentage(prunedLinks, totalLinks))
}

// updateScore takes a potential pair of contigs and update scores
func updateScore(at, bt string, score int, ctgToAlleleGroup map[string][]int, scores map[CtgAlleleGroupPair]int) {
	if gg, ok := ctgToAlleleGroup[bt]; ok {
		// Update through all alleleGroups that contig b sits in
		for _, bg := range gg {
			pair := CtgAlleleGroupPair{at, bg}
			// TODO: the score should ideally be 'normalized' score, since contig size can affect size, and as a
			//       result, the "best match" in absolute score may be wrong
			if sc, ok := scores[pair]; ok {
				if sc < score {
					scores[pair] = score
				}
			} else {
				scores[pair] = score
			}
		}
	}
}

// getScore takes a pair of contigs and get the maximum score to the allele group
func getScore(at, bt string, ctgToAlleleGroup map[string][]int, scores map[CtgAlleleGroupPair]int) int {
	if gg, ok := ctgToAlleleGroup[bt]; ok {
		maxScore := -1
		for _, bg := range gg {
			maxScoreAlleleGroup, _ := scores[CtgAlleleGroupPair{at, bg}]
			if maxScore < maxScoreAlleleGroup {
				maxScore = maxScoreAlleleGroup
			}
		}
		return maxScore // maximum score among all allele group matchings
	}
	return -1
}

// parseAllelesFile() routes parser to either parseAssociationLog() or
// parseAllelesTable(), based on the first row of the file
func parseAllelesFile(filename string) ([]AlleleGroup, error) {
	fh, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	reader := bufio.NewReader(fh)
	row, err := reader.ReadString('\n')
	if err != nil {
		return nil, err
	}
	if strings.Contains(row, "->") {
		return parseAssociationLog(filename)
	}
	err = fh.Close()
	if err != nil {
		return nil, err
	}
	return parseAllelesTable(filename)
}

// parseAssociationLog imports contig allelic relationship from purge-haplotigs
// File has the following format:
// tig00030660,PRIMARY -> tig00003333,HAPLOTIG
//                     -> tig00038686,HAPLOTIG
func parseAssociationLog(associationFile string) ([]AlleleGroup, error) {
	log.Noticef("Parse association log `%s`", associationFile)
	fh, err := os.Open(associationFile)
	if err != nil {
		return nil, err
	}

	reader := bufio.NewReader(fh)
	var data []AlleleGroup
	var primary string
	var haplotig string
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if err == io.EOF {
			break
		}
		if row == "" {
			continue
		}
		if err != nil {
			return nil, err
		}
		words := strings.Split(row, " ")
		if len(words) == 3 {
			primary = strings.Split(words[0], ",")[0]
			haplotig = strings.Split(words[2], ",")[0]
		} else if len(words) == 2 {
			haplotig = strings.Split(words[1], ",")[0]
		} else {
			return nil, fmt.Errorf("malformed line: %s, expecting 2 or 3 words", row)
		}
		alleleGroup := AlleleGroup{primary, haplotig}
		data = append(data, alleleGroup)
	}

	err = fh.Close()
	return data, err
}

// parseAllelesTable imports the contig allelic relationship
// File is a tab-separated file that looks like the following:
// Chr10   18902   tig00030660     tig00003333
// Chr10   35071   tig00038687     tig00038686     tig00065419
func parseAllelesTable(allelesFile string) ([]AlleleGroup, error) {
	log.Noticef("Parse alleles table `%s`", allelesFile)
	fh, err := os.Open(allelesFile)
	if err != nil {
		return nil, err
	}

	reader := bufio.NewReader(fh)

	var data []AlleleGroup
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		words := strings.Split(row, "\t")
		if len(words) <= 3 { // Must have at least 4 fields, i.e. 1 pair
			continue
		}
		alleleGroup := make(AlleleGroup, len(words)-2)
		copy(alleleGroup, words[2:])
		data = append(data, alleleGroup)
	}

	err = fh.Close()
	return data, err
}

// writePairsFile simply writes pruned contig pairs to file
func writePairsFile(pairsFile string, edges []ContigPair) error {
	f, _ := os.Create(pairsFile)
	w := bufio.NewWriter(f)
	_, err := fmt.Fprintf(w, PairsFileHeader)
	if err != nil {
		return err
	}

	for _, c := range edges {
		_, err = fmt.Fprintln(w, c)
		if err != nil {
			return err
		}
	}
	err = w.Flush()
	if err != nil {
		return err
	}
	log.Noticef("Pruned contig pair analyses written to `%s`", pairsFile)
	return f.Close()
}
