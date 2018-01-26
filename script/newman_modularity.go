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
	"fmt"
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
	a, b   int
	weight int
}

// Make2DSlice allocates a 2D matrix with shape (m, n)
func Make2DSlice(m, n int) [][]int {
	P := make([][]int, m)
	for i := 0; i < m; i++ {
		P[i] = make([]int, n)
	}
	return P
}

// ParseGraph imports a graph in its "edge list" form. Returns an adjacency matrix.
// a b weight
func ParseGraph(filename string) [][]int {
	file, _ := os.Open(filename)
	reader := bufio.NewReader(file)
	var edges []Edge
	var nodes []string
	nodesIdx := make(map[string]int)

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
		ai, aok := nodesIdx[a]
		if !aok {
			ai = len(nodes)
			nodesIdx[a] = ai
			nodes = append(nodes, a)
		}
		bi, bok := nodesIdx[b]
		if !bok {
			bi = len(nodes)
			nodesIdx[b] = bi
			nodes = append(nodes, b)
		}

		if len(words) == 2 {
			edges = append(edges, Edge{ai, bi, 1})
		} else {
			weight, _ := strconv.Atoi(words[2])
			edges = append(edges, Edge{ai, bi, weight})
		}
	}

	m := len(edges)
	n := len(nodes)
	log.Noticef("Graph contains %d nodes and %d edges", m, n)

	A := Make2DSlice(n, n)
	for _, e := range edges {
		A[e.a][e.b] = e.weight
		A[e.b][e.a] = e.weight
	}
	return A
}

func main() {
	logging.SetBackend(BackendFormatter)
	A := ParseGraph("karate/out.ucidata-zachary")
	fmt.Println(A)
}
