/*
 * Filename: /Users/bao/code/allhic/anchor.go
 * Path: /Users/bao/code/allhic
 * Created Date: Monday, June 4th 2018, 9:26:26 pm
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
	"sort"

	"github.com/biogo/hts/bam"
)

// Anchorer runs the merging algorithm
type Anchorer struct {
	Bamfile      string
	contigs      []Contig
	nameToContig map[string]int
	links        [][]Link
}

// Contig stores the name and length of each contig
type Contig struct {
	id     int
	name   string
	length int
}

// Link contains a specific inter-contig link
type Link struct {
	a, b       int // Contig ids
	apos, bpos int // Positions
}

// Run kicks off the merging algorithm
func (r *Anchorer) Run() {
	r.ExtractInterContigLinks()
	log.Notice("Success")
}

// ExtractInterContigLinks extracts links from the Bamfile
func (r *Anchorer) ExtractInterContigLinks() {
	fh, _ := os.Open(r.Bamfile)
	prefix := RemoveExt(r.Bamfile)
	disfile := prefix + ".dis"
	idsfile := prefix + ".ids"

	log.Noticef("Parse bamfile `%s`", r.Bamfile)
	br, _ := bam.NewReader(fh, 0)
	defer br.Close()

	fdis, _ := os.Create(disfile)
	wdis := bufio.NewWriter(fdis)
	fids, _ := os.Create(idsfile)
	wids := bufio.NewWriter(fids)

	r.nameToContig = make(map[string]int)
	refs := br.Header().Refs()
	for _, ref := range refs {
		contig := Contig{
			id:     len(r.contigs),
			name:   ref.Name(),
			length: ref.Len(),
		}
		r.contigs = append(r.contigs, contig)
		r.nameToContig[contig.name] = contig.id
		fmt.Fprintf(wids, "%s\t%d\n", ref.Name(), ref.Len())
	}
	wids.Flush()
	log.Noticef("Extracted %d contigs to `%s`", len(r.contigs), idsfile)

	// Import links into pairs of contigs
	total := 0
	intraLinks := make(map[string][]int)
	r.links = make([][]Link, len(r.contigs))
	for {
		rec, err := br.Read()
		if err != nil {
			if err != io.EOF {
				log.Error(err)
			}
			break
		}

		at, bt := rec.Ref.Name(), rec.MateRef.Name()
		a, b := r.nameToContig[at], r.nameToContig[bt]
		apos, bpos := rec.Pos, rec.MatePos

		// An intra-contig link
		if a == b {
			if link := abs(apos - bpos); link >= MinLinkDist {
				intraLinks[at] = append(intraLinks[at], link)
			}
			continue
		}

		//         read1                                               read2
		//     ---a-- X|----- dist = a2 ----|         |--- dist = b ---|X ------ b2 ------
		//     ==============================         ====================================
		//             C1 (length L1)       |----D----|         C2 (length L2)
		// An inter-contig link

		r.links[a] = append(r.links[a], Link{
			a: a, b: b, apos: apos, bpos: bpos,
		})
	}

	for _, links := range r.links {
		sort.Slice(links, func(i, j int) bool {
			return links[i].apos < links[j].apos
		})
	}

	for i, links := range r.links {
		fmt.Println(i, ":", links)
	}

	// Write intra-links to .dis file
	for contig, links := range intraLinks {
		links = unique(links)
		total += len(links)
		fmt.Fprintf(wdis, "%s\t%s\n", contig, arrayToString(links, ","))
	}
	wdis.Flush()
	log.Noticef("Extracted %d intra-contig link groups to `%s` (total = %d)",
		len(intraLinks), disfile, total)
}
