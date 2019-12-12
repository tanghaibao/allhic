/**
 * Filename: /Users/htang/code/allhic/main.go
 * Path: /Users/htang/code/allhic
 * Created Date: Wednesday, January 3rd 2018, 11:21:45 am
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"fmt"
	"github.com/spf13/cobra"
	"path"
	"sort"
	"strconv"
	"strings"
)

// Logo banner (Varsity style):
// http://patorjk.com/software/taag/#p=testall&f=3D-ASCII&t=ALLHIC
var rootCmd = &cobra.Command{
	Use:   "allhic",
	Short: "Genome scaffolding based on Hi-C data",
	Long: `
     _       _____     _____     ____  ____  _____   ______
    / \     |_   _|   |_   _|   |_   ||   _||_   _|.' ___  |
   / _ \      | |       | |       | |__| |    | | / .'   \_|
  / ___ \     | |   _   | |   _   |  __  |    | | | |
_/ /   \ \_  _| |__/ | _| |__/ | _| |  | |_  _| |_\ ` + "`" + `.___.'\
|____| |____||________||________||____||____||_____|` + "`" + `.____ .'

`,
	Version: Version,
}

// Execute executes the root command.
func Execute() error {
	return rootCmd.Execute()
}

// banner prints the separate steps
func banner(message string) {
	message = "* " + message + " *"
	log.Noticef(strings.Repeat("*", len(message)))
	log.Noticef(message)
	log.Noticef(strings.Repeat("*", len(message)))
}

// init adds all the sub-commands
func init() {
	var RE string
	var minLinks int
	extractCmd := &cobra.Command{
		Use:   "extract bamfile fastafile",
		Short: "Extract Hi-C link size distribution",
		Long: `
Extract function:
Given a bamfile, the goal of the extract step is to calculate an empirical
distribution of Hi-C link size based on intra-contig links. The Extract function
also prepares for the latter steps of ALLHiC.
`,
		Args: cobra.ExactArgs(2),
		Run: func(cmd *cobra.Command, args []string) {
			bamfile := args[0]
			fastafile := args[1]
			p := Extracter{Bamfile: bamfile, Fastafile: fastafile, RE: RE, MinLinks: minLinks}
			p.Run()
		},
	}
	extractCmd.Flags().StringVarP(&RE, "RE", "", DefaultRE, "Restriction site pattern, use comma to separate multiple patterns (N is considered as [ACGT]), e.g. 'GATCGATC,GANTGATC,GANTANTC,GATCANTC'")
	extractCmd.Flags().IntVarP(&minLinks, "minLinks", "", MinLinks, "Minimum number of links for contig pair")

	allelesCmd := &cobra.Command{
		Use:   "alleles genome.paf genome.counts_RE.txt",
		Short: "Build alleles.table for `prune`",
		Long: `
Alleles function:
Given a paf file, we could identify and classify the allelic contigs to be used
for "allhic prune". We recommend the following parameters to build the paf file:

$ minimap2 -DP -k19 -w19 -m200 -t32 genome.fasta genome.fasta > genome.paf

The PAF file contains all self-alignments, which is the basis for classification.
ALLHiC generates "alleles.table", which can then be used for later steps.
`,
		Args: cobra.ExactArgs(2),
		Run: func(cmd *cobra.Command, args []string) {
			pafFile := args[0]
			reFile := args[1]
			p := Alleler{PafFile: pafFile, ReFile: reFile}
			p.Run()
		},
	}

	pruneCmd := &cobra.Command{
		Use:   "prune alleles.table pairs.txt",
		Short: "Prune allelic, cross-allelic and weak links",
		Long: `
Prune function:
Given contig pairing, the goal of the prune step is to remove all inter-allelic
links, then it is possible to reconstruct allele-separated assemblies in the following
steps. The alleles.table file contains tab-separated columns containing contigs that
are considered "allelic". For example:

Chr10   18902   tig00030660     tig00003333
Chr10   35071   tig00038687     tig00038686     tig00065419

These lines means that at certain locations in the genome (typically based on synteny
to a closely-related genome), several contigs are considered allelic. The Hi-C links
between these contigs are removed, in addition, other contigs linking to these contigs
simultaneously would only consider one single-best edge. The "pairs.txt" file is the output
of the "extract" command.

Alternatively we can use output from purge-haplotigs as source of alleles.table, the format
looks like:

tig00030660,PRIMARY -> tig00003333,HAPLOTIG
                    -> tig00038686,HAPLOTIG
`,
		Args: cobra.ExactArgs(2),
		Run: func(cmd *cobra.Command, args []string) {
			allelesFile := args[0]
			pairsFile := args[1]
			p := Pruner{AllelesFile: allelesFile, PairsFile: pairsFile}
			p.Run()
		},
	}

	var minREs, maxLinkDensity, nonInformativeRatio int
	partitionCmd := &cobra.Command{
		Use:   "partition counts_RE.txt pairs.txt k",
		Short: "Separate contigs into k groups",
		Long: `
Partition function:
Given a target k, number of partitions, the goal of the partitioning is to
separate all the contigs into separate clusters. As with all clustering
algorithm, there is an optimization goal here. The LACHESIS algorithm is
a hierarchical clustering algorithm using average links. The two input files
can be generated with the "extract" sub-command.
`,
		Args: cobra.ExactArgs(3),
		Run: func(cmd *cobra.Command, args []string) {
			contigsfile := args[0]
			pairsFile := args[1]
			k, _ := strconv.Atoi(args[2])
			p := Partitioner{Contigsfile: contigsfile, PairsFile: pairsFile, K: k,
				MinREs: minREs, MaxLinkDensity: maxLinkDensity,
				NonInformativeRatio: nonInformativeRatio}
			p.Run()
		},
	}
	partitionCmd.Flags().IntVarP(&minREs, "minREs", "", MinREs, "Minimum number of RE sites in a contig to be clustered (CLUSTER_MIN_RE_SITES in LACHESIS)")
	partitionCmd.Flags().IntVarP(&maxLinkDensity, "maxLinkDensity", "", MaxLinkDensity, "Density threshold before marking contig as repetive (CLUSTER_MAX_LINK_DENSITY in LACHESIS)")
	partitionCmd.Flags().IntVarP(&nonInformativeRatio, "nonInformativeRatio", "", NonInformativeRatio, "cutoff for recovering skipped contigs back into the clusters (CLUSTER_NONINFORMATIVE_RATIO in LACHESIS)")

	var skipGA, resume bool
	var seed int64
	var npop, ngen int
	var mutpb float64
	optimizeCmd := &cobra.Command{
		Use:   "optimize counts_RE.txt clmfile",
		Short: "Order-and-orient tigs in a group",
		Long: `
Optimize function:
Given a set of Hi-C contacts between contigs, as specified in the
clmfile, reconstruct the highest scoring ordering and orientations
for these contigs. Optimize run on a specific partition in "clusters.txt"
as generated by "partition" sub-command, with the group_number matching the
order appearing in "clusters.txt". Typically, if there are k clusters, we
can start k separate "optimize" commands for parallelism (for example,
on a cluster).
`,
		Args: cobra.ExactArgs(2),
		Run: func(cmd *cobra.Command, args []string) {
			refile := args[0]
			clmfile := args[1]
			p := Optimizer{REfile: refile, Clmfile: clmfile,
				RunGA: !skipGA, Resume: resume,
				Seed: seed, NPop: npop, NGen: ngen, MutProb: mutpb}
			p.Run()
		},
	}
	optimizeCmd.Flags().BoolVarP(&skipGA, "skipGA", "", false, "Skip GA step")
	optimizeCmd.Flags().BoolVarP(&resume, "resume", "", false, "Resume from existing tour file")
	optimizeCmd.Flags().Int64VarP(&seed, "seed", "", Seed, "Random seed")
	optimizeCmd.Flags().IntVarP(&npop, "npop", "", Npop, "Population size")
	optimizeCmd.Flags().IntVarP(&ngen, "ngen", "", Ngen, "Number of generations for convergence")
	optimizeCmd.Flags().Float64VarP(&mutpb, "mutapb", "", MutaProb, "Mutation prob in GA")

	buildCmd := &cobra.Command{
		Use:   "build tourfile1 tourfile2 ... contigs.fasta asm.chr.fasta",
		Short: "Build genome release",
		Long: `
Build function:
Convert the tourfile into the standard AGP file, which is then converted
into a FASTA genome release.
`,
		Args: cobra.MinimumNArgs(3),
		Run: func(cmd *cobra.Command, args []string) {

			tourfiles := []string{}
			for i := 0; i < len(args)-2; i++ {
				tourfiles = append(tourfiles, args[i])
			}
			sort.Strings(tourfiles)

			fastafile := args[len(args)-2]
			outfastafile := args[len(args)-1]
			p := Builder{Tourfiles: tourfiles,
				Fastafile:    fastafile,
				OutFastafile: outfastafile}
			p.Run()
		},
	}

	plotCmd := &cobra.Command{
		Use:   "plot bamfile tourfile",
		Short: "Extract matrix of link counts and plot heatmap",
		Long: `
Anchor function:
Given a bamfile, we extract matrix of link counts and plot heatmap.
`,
		Args: cobra.ExactArgs(2),
		Run: func(cmd *cobra.Command, args []string) {
			bamfile := args[0]
			tourfile := args[1]
			p := Plotter{
				Anchor: &Anchorer{Bamfile: bamfile, Tourfile: tourfile},
			}
			p.Run()
		},
	}

	assessCmd := &cobra.Command{
		Use:   "assess bamfile bedfile chr1",
		Short: "Assess the orientations of contigs",
		Long: `
Assess function:
Compute the posterior probability of contig orientations after scaffolding
as a quality assessment step.
`,
		Args: cobra.ExactArgs(3),
		Run: func(cmd *cobra.Command, args []string) {
			bamfile := args[0]
			bedfile := args[1]
			seqid := args[2]
			p := Assesser{Bamfile: bamfile, Bedfile: bedfile, Seqid: seqid}
			p.Run()
		},
	}

	pipelineCmd := &cobra.Command{
		Use:   "pipeline bamfile fastafile k",
		Short: "Run extract-partition-optimize-build steps sequentially",
		Long: `
Pipeline:
A convenience driver function. Chain the following steps sequentially.

- extract
- partion
- optimize
- build
`,
		Args: cobra.ExactArgs(3),
		Run: func(cmd *cobra.Command, args []string) {
			bamfile := args[0]
			fastafile := args[1]
			k, _ := strconv.Atoi(args[2])

			// Extract the contig pairs, count RE sites
			banner(fmt.Sprintf("Extractor started (RE = %s)", RE))
			extractor := Extracter{Bamfile: bamfile, Fastafile: fastafile, RE: RE}
			extractor.Run()

			// Partition into k groups
			banner(fmt.Sprintf("Partition into %d groups", k))
			partitioner := Partitioner{Contigsfile: extractor.OutContigsfile,
				PairsFile: extractor.OutPairsfile, K: k}
			partitioner.Run()

			// Optimize the k groups separately
			tourfiles := []string{}
			for i, refile := range partitioner.OutREfiles {
				banner(fmt.Sprintf("Optimize group %d", i))
				optimizer := Optimizer{REfile: refile,
					Clmfile: extractor.OutClmfile,
					RunGA:   !skipGA, Resume: resume,
					Seed: seed, NPop: npop, NGen: ngen, MutProb: mutpb}
				optimizer.Run()
				tourfiles = append(tourfiles, optimizer.OutTourFile)
			}

			// Run the final build
			banner("Build started (AGP and FASTA)")
			outfastafile := path.Join(path.Dir(tourfiles[0]),
				fmt.Sprintf("asm-g%d.chr.fasta", k))
			builder := Builder{Tourfiles: tourfiles,
				Fastafile:    fastafile,
				OutFastafile: outfastafile}
			builder.Run()
		},
	}
	pipelineCmd.Flags().StringVarP(&RE, "RE", "", DefaultRE, "Restriction site pattern, use comma to separate multiple patterns (N is considered as [ACGT]), e.g. 'GATCGATC,GANTGATC,GANTANTC,GATCANTC'")
	pipelineCmd.Flags().IntVarP(&minLinks, "minLinks", "", MinLinks, "Minimum number of links for contig pair")

	pipelineCmd.Flags().IntVarP(&minREs, "minREs", "", MinREs, "Minimum number of RE sites in a contig to be clustered (CLUSTER_MIN_RE_SITES in LACHESIS)")
	pipelineCmd.Flags().IntVarP(&maxLinkDensity, "maxLinkDensity", "", MaxLinkDensity, "Density threshold before marking contig as repetive (CLUSTER_MAX_LINK_DENSITY in LACHESIS)")
	pipelineCmd.Flags().IntVarP(&nonInformativeRatio, "nonInformativeRatio", "", NonInformativeRatio, "cutoff for recovering skipped contigs back into the clusters (CLUSTER_NONINFORMATIVE_RATIO in LACHESIS)")

	pipelineCmd.Flags().BoolVarP(&skipGA, "skipGA", "", false, "Skip GA step")
	pipelineCmd.Flags().BoolVarP(&resume, "resume", "", false, "Resume from existing tour file")
	pipelineCmd.Flags().Int64VarP(&seed, "seed", "", Seed, "Random seed")
	pipelineCmd.Flags().IntVarP(&npop, "npop", "", Npop, "Population size")
	pipelineCmd.Flags().IntVarP(&ngen, "ngen", "", Ngen, "Number of generations for convergence")
	pipelineCmd.Flags().Float64VarP(&mutpb, "mutapb", "", MutaProb, "Mutation prob in GA")

	rootCmd.AddCommand(extractCmd, allelesCmd, pruneCmd, partitionCmd, optimizeCmd, buildCmd, plotCmd, assessCmd, pipelineCmd)
}
