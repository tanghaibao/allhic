/**
 * Filename: /Users/htang/code/allhic/main.go
 * Path: /Users/htang/code/allhic
 * Created Date: Wednesday, January 3rd 2018, 11:21:45 am
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

package main

import (
	"fmt"
	"os"
	"path"
	"sort"
	"strconv"
	"strings"
	"time"

	logging "github.com/op/go-logging"
	"github.com/tanghaibao/allhic"
	"github.com/urfave/cli"
)

var log = logging.MustGetLogger("main")

// init customizes how cli layout the command interface
// Logo banner (Varsity style):
// http://patorjk.com/software/taag/#p=testall&f=3D-ASCII&t=ALLHIC
func init() {
	cli.AppHelpTemplate = `
     _       _____     _____     ____  ____  _____   ______
    / \     |_   _|   |_   _|   |_   ||   _||_   _|.' ___  |
   / _ \      | |       | |       | |__| |    | | / .'   \_|
  / ___ \     | |   _   | |   _   |  __  |    | | | |
_/ /   \ \_  _| |__/ | _| |__/ | _| |  | |_  _| |_\ ` + "`" + `.___.'\
|____| |____||________||________||____||____||_____|` + "`" + `.____ .'

` + cli.AppHelpTemplate
}

// banner prints the separate steps
func banner(message string) {
	message = "* " + message + " *"
	log.Noticef(strings.Repeat("*", len(message)))
	log.Noticef(message)
	log.Noticef(strings.Repeat("*", len(message)))
}

// main is the entrypoint for the entire program, routes to commands
func main() {
	logging.SetBackend(allhic.BackendFormatter)

	app := cli.NewApp()
	app.Compiled = time.Now()
	app.Copyright = "(c) Haibao Tang, Xingtan Zhang 2017-2018"
	app.Name = "ALLHIC"
	app.Usage = "Genome scaffolding based on Hi-C data"
	app.Version = allhic.Version

	extractFlags := []cli.Flag{
		cli.StringFlag{
			Name:  "RE",
			Usage: "Restriction site pattern",
			Value: "GATC",
		},
	}

	partitionFlags := []cli.Flag{
		cli.IntFlag{
			Name:  "minREs",
			Usage: "Minimum number of RE sites in a contig to be clustered (CLUSTER_MIN_RE_SITES in LACHESIS)",
			Value: allhic.MinREs,
		},
		cli.IntFlag{
			Name:  "maxLinkDensity",
			Usage: "Density threshold before marking contig as repetive (CLUSTER_MAX_LINK_DENSITY in LACHESIS)",
			Value: allhic.MaxLinkDensity,
		},
		cli.IntFlag{
			Name:  "nonInformativeRatio",
			Usage: "cutoff for recovering skipped contigs back into the clusters (CLUSTER_NONINFORMATIVE_RATIO in LACHESIS)",
			Value: allhic.NonInformativeRatio,
		},
	}

	optimizeFlags := []cli.Flag{
		cli.BoolFlag{
			Name:  "skipGA",
			Usage: "Skip GA step",
		},
		cli.BoolFlag{
			Name:  "resume",
			Usage: "Resume from existing tour file",
		},
		cli.Int64Flag{
			Name:  "seed",
			Usage: "Random seed",
			Value: 42,
		},
		cli.IntFlag{
			Name:  "npop",
			Usage: "Population size",
			Value: 100,
		},
		cli.IntFlag{
			Name:  "ngen",
			Usage: "Number of generations for convergence",
			Value: 5000,
		},
		cli.Float64Flag{
			Name:  "mutpb",
			Usage: "Mutation prob in GA",
			Value: .2,
		},
	}

	app.Commands = []cli.Command{
		{
			Name:  "extract",
			Usage: "Extract Hi-C link size distribution",
			UsageText: `
	allhic extract bamfile fastafile [options]

Extract function:
Given a bamfile, the goal of the extract step is to calculate an empirical
distribution of Hi-C link size based on intra-contig links. The Extract function
also prepares for the latter steps of ALLHIC.
`,
			Flags: extractFlags,
			Action: func(c *cli.Context) error {
				if c.NArg() < 2 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify bamfile and fastafile", 1)
				}

				bamfile := c.Args().Get(0)
				fastafile := c.Args().Get(1)
				RE := c.String("RE")

				p := allhic.Extracter{Bamfile: bamfile, Fastafile: fastafile, RE: RE}
				p.Run()
				return nil
			},
		},
		{
			Name:  "alleles",
			Usage: "Build alleles.table for `prune`",
			UsageText: `
	allhic alleles genome.paf genome.counts_RE.txt

Alleles function:
Given a paf file, we could identify and classify the allelic contigs to be used
for "allhic prune". We recommend the following parameters to build the paf file:

$ minimap2 -DP -k19 -w19 -m200 -t32 genome.fasta genome.fasta > genome.paf

The PAF file contains all self-alignments, which is the basis for classification.
ALLHiC generates "alleles.table", which can then be used for later steps.
`,
			Action: func(c *cli.Context) error {
				if c.NArg() < 2 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specificify paf file", 1)
				}

				pafFile := c.Args().Get(0)
				reFile := c.Args().Get(1)
				p := allhic.Alleler{PafFile: pafFile, REFile: reFile}
				p.Run()
				return nil
			},
		},
		{
			Name:  "prune",
			Usage: "Prune allelic, cross-allelic and weak links",
			UsageText: `
	allhic prune alleles.table pairs.txt [options]

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
`,
			Action: func(c *cli.Context) error {
				if c.NArg() < 2 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify alleles.table and pairs.txt", 1)
				}

				allelesFile := c.Args().Get(0)
				pairsFile := c.Args().Get(1)
				p := allhic.Pruner{AllelesFile: allelesFile, PairsFile: pairsFile}
				p.Run()
				return nil
			},
		},
		{
			Name:  "partition",
			Usage: "Separate contigs into k groups",
			UsageText: `
	allhic partition counts_RE.txt pairs.txt k [options]

Partition function:
Given a target k, number of partitions, the goal of the partitioning is to
separate all the contigs into separate clusters. As with all clustering
algorithm, there is an optimization goal here. The LACHESIS algorithm is
a hierarchical clustering algorithm using average links. The two input files
can be generated with the "extract" sub-command.
`,
			Flags: partitionFlags,
			Action: func(c *cli.Context) error {
				if c.NArg() < 3 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify countsFile, pairsFile and value k", 1)
				}

				contigsfile := c.Args().Get(0)
				pairsFile := c.Args().Get(1)
				k, _ := strconv.Atoi(c.Args().Get(2))
				minREs := c.Int("minREs")
				maxLinkDensity := c.Int("maxLinkDensity")
				nonInformativeRatio := c.Int("nonInformativeRatio")
				p := allhic.Partitioner{Contigsfile: contigsfile, PairsFile: pairsFile, K: k,
					MinREs: minREs, MaxLinkDensity: maxLinkDensity,
					NonInformativeRatio: nonInformativeRatio}
				p.Run()
				return nil
			},
		},
		{
			Name:  "optimize",
			Usage: "Order-and-orient tigs in a group",
			UsageText: `
	allhic optimize counts_RE.txt clmfile [options]

Optimize function:
Given a set of Hi-C contacts between contigs, as specified in the
clmfile, reconstruct the highest scoring ordering and orientations
for these contigs. Optimize run on a specific partition in "clusters.txt"
as generated by "partition" sub-command, with the group_number matching the
order appearing in "clusters.txt". Typically, if there are k clusters, we
can start k separate "optimize" commands for parallelism (for example,
on a cluster).
`,
			Flags: optimizeFlags,
			Action: func(c *cli.Context) error {
				if c.NArg() < 2 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify clmfile", 1)
				}

				refile := c.Args().Get(0)
				clmfile := c.Args().Get(1)
				runGA := !c.Bool("skipGA")
				resume := c.Bool("resume")
				seed := c.Int64("seed")
				npop := c.Int("npop")
				ngen := c.Int("ngen")
				mutpb := c.Float64("mutpb")
				p := allhic.Optimizer{REfile: refile, Clmfile: clmfile,
					RunGA: runGA, Resume: resume,
					Seed: seed, NPop: npop, NGen: ngen, MutProb: mutpb}
				p.Run()
				return nil
			},
		},
		{
			Name:  "build",
			Usage: "Build genome release",
			UsageText: `
	allhic build tourfile1 tourfile2 ... contigs.fasta asm.chr.fasta  [options]

Build function:
Convert the tourfile into the standard AGP file, which is then converted
into a FASTA genome release.
`,
			Action: func(c *cli.Context) error {
				if c.NArg() < 3 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify tourfile and fastafile", 1)
				}

				tourfiles := []string{}
				for i := 0; i < c.NArg()-2; i++ {
					tourfiles = append(tourfiles, c.Args().Get(i))
				}
				sort.Strings(tourfiles)

				fastafile := c.Args().Get(c.NArg() - 2)
				outfastafile := c.Args().Get(c.NArg() - 1)

				p := allhic.Builder{Tourfiles: tourfiles,
					Fastafile:    fastafile,
					OutFastafile: outfastafile}
				p.Run()
				return nil
			},
		},
		{
			Name:  "plot",
			Usage: "Extract matrix of link counts and plot heatmap",
			UsageText: `
	allhic anchor bamfile tourfile [options]

Anchor function:
Given a bamfile, we extract matrix of link counts and plot heatmap.
`,
			Action: func(c *cli.Context) error {
				if c.NArg() < 2 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify bamfile", 1)
				}

				bamfile := c.Args().Get(0)
				tourfile := c.Args().Get(1)
				p := allhic.Plotter{
					Anchor: &allhic.Anchorer{Bamfile: bamfile, Tourfile: tourfile},
				}
				p.Run()
				return nil
			},
		},
		{
			Name:  "assess",
			Usage: "Assess the orientations of contigs",
			UsageText: `
	allhic assess bamfile bedfile chr1

Assess function:
Compute the posterior probability of contig orientations after scaffolding
as a quality assessment step.
`,
			Action: func(c *cli.Context) error {
				if c.NArg() < 3 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify bamfile", 1)
				}

				bamfile := c.Args().Get(0)
				bedfile := c.Args().Get(1)
				seqid := c.Args().Get(2)
				p := allhic.Assesser{Bamfile: bamfile, Bedfile: bedfile, Seqid: seqid}
				p.Run()
				return nil
			},
		},
		{
			Name:  "pipeline",
			Usage: "Run extract-partition-optimize-build steps sequentially",
			UsageText: `
	allhic pipeline bamfile fastafile k [options]

Pipeline:
A convenience driver function. Chain the following steps sequentially.

- extract
- partion
- optimize
- build
`,
			Flags: append(extractFlags, optimizeFlags...),
			Action: func(c *cli.Context) error {
				if c.NArg() < 3 {
					cli.ShowSubcommandHelp(c)
					return cli.NewExitError("Must specify pairsFile, clmfile and bamfile", 1)
				}

				bamfile := c.Args().Get(0)
				fastafile := c.Args().Get(1)
				k, _ := strconv.Atoi(c.Args().Get(2))
				RE := c.String("RE")
				runGA := !c.Bool("skipGA")
				resume := c.Bool("resume")
				seed := c.Int64("seed")
				npop := c.Int("npop")
				ngen := c.Int("ngen")
				mutpb := c.Float64("mutpb")

				// Extract the contig pairs, count RE sites
				banner(fmt.Sprintf("Extractor started (RE = %s)", RE))
				extractor := allhic.Extracter{Bamfile: bamfile, Fastafile: fastafile, RE: RE}
				extractor.Run()

				// Partition into k groups
				banner(fmt.Sprintf("Partition into %d groups", k))
				partitioner := allhic.Partitioner{Contigsfile: extractor.OutContigsfile,
					PairsFile: extractor.OutPairsfile, K: k}
				partitioner.Run()

				// Optimize the k groups separately
				tourfiles := []string{}
				for i, refile := range partitioner.OutREfiles {
					banner(fmt.Sprintf("Optimize group %d", i))
					optimizer := allhic.Optimizer{REfile: refile,
						Clmfile: extractor.OutClmfile,
						RunGA:   runGA, Resume: resume,
						Seed: seed, NPop: npop, NGen: ngen, MutProb: mutpb}
					optimizer.Run()
					tourfiles = append(tourfiles, optimizer.OutTourFile)
				}

				// Run the final build
				banner("Build started (AGP and FASTA)")
				outfastafile := path.Join(path.Dir(tourfiles[0]),
					fmt.Sprintf("asm-g%d.chr.fasta", k))
				builder := allhic.Builder{Tourfiles: tourfiles,
					Fastafile:    fastafile,
					OutFastafile: outfastafile}
				builder.Run()

				return nil
			},
		},
	}

	app.Run(os.Args)
}
