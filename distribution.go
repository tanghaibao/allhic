/**
 * Filename: /Users/bao/code/allhic/distribution.go
 * Path: /Users/bao/code/allhic
 * Created Date: Wednesday, March 7th 2018, 1:56:45 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"
)

// Distribution processes the distribution step
type Distribution struct {
	Bamfile     string
	minLinkDist int
	maxLinkDist int
	binStarts   []int
	links       []int
}

// Run calls the distribution steps
func (r *Distribution) Run() {
	file, _ := os.Open(r.Bamfile)
	log.Noticef("Parse dist file `%s`", r.Bamfile)
	reader := bufio.NewReader(file)
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		words := strings.Split(row, "\t")
		for _, link := range strings.Split(words[1], ",") {
			ll, _ := strconv.Atoi(link)
			if ll > 0 {
				r.links = append(r.links, ll)
			}
		}
	}
	log.Noticef("Imported %d intra-contig links", len(r.links))

	r.Makebins()
}

// Makebins makes geometric bins and count links that fall in each bin
// This borrows the method form LACHESIS
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeDistribution.cc
func (r *Distribution) Makebins() {
	r.minLinkDist, r.maxLinkDist = MaxInt, MinInt
	for _, link := range r.links {
		if link > r.maxLinkDist {
			r.maxLinkDist = link
		}
		if link < r.minLinkDist {
			r.minLinkDist = link
		}
	}
	linkRange := math.Log2(float64(r.maxLinkDist) / float64(r.minLinkDist))
	nBins := int(math.Ceil(linkRange * 16))
	log.Noticef("Min link: %d, Max link: %d, # Bins: %d",
		r.minLinkDist, r.maxLinkDist, nBins)

	for i := 0; 16*i <= nBins; i++ {
		jpower := 1.0
		for j := 0; j < 16 && 16*i+j <= nBins; j++ {
			binStart := float64(r.minLinkDist<<uint(i)) * jpower
			r.binStarts = append(r.binStarts, int(binStart))
			jpower *= GeometricBinSize
		}
	}
	fmt.Println(r.binStarts)
}
