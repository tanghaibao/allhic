/**
 * Filename: /Users/bao/code/allhic/allhic/clm.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Monday, January 1st 2018, 5:57:00 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"bufio"
	"fmt"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/op/go-logging"
)

var log = logging.MustGetLogger("example")
var format = logging.MustStringFormatter(
	`%{color}%{time:15:04:05.000} %{shortfunc} ▶ %{level:.4s} %{color:reset} %{message}`,
)
var backend = logging.NewLogBackend(os.Stderr, "", 0)
var backendFormatter = logging.NewBackendFormatter(backend, format)

// RemoveExt returns the substring minus the extension
func RemoveExt(filename string) string {
	return strings.TrimSuffix(filename, path.Ext(filename))
}

// CLMFile has the following format:
//
// tig00046211+ tig00063795+       1       53173
// tig00046211+ tig00063795-       1       116050
// tig00046211- tig00063795+       1       71155
// tig00046211- tig00063795-       1       134032
// tig00030676+ tig00077819+       7       136407 87625 87625 106905 102218 169660 169660
// tig00030676+ tig00077819-       7       126178 152952 152952 35680 118923 98367 98367
// tig00030676- tig00077819+       7       118651 91877 91877 209149 125906 146462 146462
// tig00030676- tig00077819-       7       108422 157204 157204 137924 142611 75169 75169
type CLMFile struct {
	Name      string
	Clmfile   string
	Idsfile   string
	tigToSize map[string]int
}

// InitCLMFile is the constructor for CLMFile
func InitCLMFile(Clmfile string) *CLMFile {
	p := new(CLMFile)
	p.Name = RemoveExt(path.Base(Clmfile))
	p.Clmfile = Clmfile
	p.Idsfile = RemoveExt(Clmfile) + ".ids"
	p.tigToSize = make(map[string]int)
	return p
}

// ParseIds parses the idsfile into data stored in CLMFile.
// IDS file has a list of contigs that need to be ordered. 'recover',
// keyword, if available in the third column, is less confident.
// tig00015093     46912
// tig00035238     46779   recover
// tig00030900     119291
func (r *CLMFile) ParseIds() {
	logging.SetBackend(backend, backendFormatter)

	file, _ := os.Open(r.Idsfile)
	log.Infof("Parse idsfile `%s`", r.Idsfile)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		words := strings.Fields(scanner.Text())
		tig := words[0]
		size, _ := strconv.Atoi(words[1])
		r.tigToSize[tig] = size
	}

	fmt.Println(r.tigToSize)
}

// ParseClm parses the clmfile into data stored in CLMFile.
func (r *CLMFile) ParseClm() {
	// var contacts map[(string, string)]int

	file, _ := os.Open(r.Clmfile)
	log.Infof("Parse clmfile `%s`", r.Clmfile)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		row := strings.TrimSpace(scanner.Text())
		words := strings.Split(row, "\t")
		abtig := strings.Split(words[0], " ")
		// atig, btig := abtig[0], abtig[1]
		// at, ao := atig[:len(atig)-1], atig[len(atig)-1]
		// bt, bo := btig[:len(btig)-1], btig[len(btig)-1]

		var dists []int
		for _, dist := range strings.Split(words[2], " ") {
			d, _ := strconv.Atoi(dist)
			dists = append(dists, d)
		}
		fmt.Println(abtig, dists)
	}
}
