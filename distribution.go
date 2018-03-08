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
	"io"
	"os"
	"strconv"
	"strings"
)

// Distribution processes the distribution step
type Distribution struct {
	Bamfile string
}

// Run calls the distribution steps
func (r *Distribution) Run() {
	file, _ := os.Open(r.Bamfile)
	log.Noticef("Parse dist file `%s`", r.Bamfile)
	reader := bufio.NewReader(file)
	var links []int
	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		words := strings.Split(row, "\t")
		for _, link := range strings.Split(words[1], ",") {
			ll, _ := strconv.Atoi(link)
			links = append(links, ll)
		}
	}
	log.Noticef("Imported %d intra-contig links", len(links))
}
