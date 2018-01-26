/**
 * Filename: /Users/htang/code/allhic/script/newman_modularity.go
 * Path: /Users/htang/code/allhic/script
 * Created Date: Thursday, January 25th 2018, 6:01:39 pm
 * Author: htang
 *
 * Copyright (c) 2018 Haibao Tang
 */

// "Prototype of the Newman modularity inference method. Implemented based on
// Newman. (2006) Modularity and community structure in networks.

package main

import (
	"bufio"
	"io"
	"os"
	"strconv"
	"strings"

	logging "github.com/op/go-logging"
)

var log = logging.MustGetLogger("allhic")
var format = logging.MustStringFormatter(
	`%{color}%{time:15:04:05} %{shortfunc} | %{level:.6s} %{color:reset} %{message}`,
)

// Backend is the default stderr output
var Backend = logging.NewLogBackend(os.Stderr, "", 0)

// BackendFormatter contains the fancy debug formatter
var BackendFormatter = logging.NewBackendFormatter(Backend, format)

// Edge represents an edge in the graph
type Edge struct {
	a, b   string
	weight int
}

// ParseGraph imports a graph in its "edge list" form. Returns an adjacency matrix.
// a b weight
func ParseGraph(filename string) {
	file, _ := os.Open(filename)
	reader := bufio.NewReader(file)
	var edges []Edge
	nodes := make(map[string]bool)

	for {
		row, err := reader.ReadString('\n')
		row = strings.TrimSpace(row)
		if row == "" && err == io.EOF {
			break
		}
		if row[0] == '%' {
			continue
		}
		words := strings.Fields(row)
		a, b := words[0], words[1]
		nodes[a] = true
		nodes[b] = true

		if len(words) == 2 {
			edges = append(edges, Edge{a, b, 1})
		} else {
			weight, _ := strconv.Atoi(words[2])
			edges = append(edges, Edge{a, b, weight})
		}
	}

	m := len(edges)
	n := len(nodes)
	log.Noticef("Graph contains %d nodes and %d edges", m, n)
}

func main() {
	logging.SetBackend(BackendFormatter)
	ParseGraph("karate/out.ucidata-zachary")
}
