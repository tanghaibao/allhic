/**
 * Filename: /Users/bao/code/allhic/prune.go
 * Path: /Users/bao/code/allhic
 * Created Date: Wednesday, February 28th 2018, 8:46:25 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

// Pruner processes the pruning step
type Pruner struct {
	Bamfile string
}

// Run calls the pruning steps
func (r *Pruner) Run() {
	log.Errorf("Prune function is still under development. Please use `./prune` for now.")
	log.Errorf("Usage: ./prune -i Allele.ctg.table -b bam.list -r draft.asm.fasta")
}
