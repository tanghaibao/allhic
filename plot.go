/*
 * Filename: /Users/bao/code/allhic/plot.go
 * Path: /Users/bao/code/allhic
 * Created Date: Saturday, July 7th 2018, 1:33:37 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import "os"

// Plotter extracts a matrix of link counts and plot a heatmp
type Plotter struct {
	Anchor *Anchorer
}

// Run starts the plotting
func (r *Plotter) Run() {
	m := r.Anchor
	m.ExtractInterContigLinks()
	m.parseTourFile(m.Tourfile)
	m.printTour(os.Stdout, "PLOTTER")
	// Serialize to disk for plotting
	m.makeContigStarts()
	m.serialize(250000, "genome.json", "data.npy")

	// printPaths(paths)
	log.Notice("Success")
}
