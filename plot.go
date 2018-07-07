/*
 * Filename: /Users/bao/code/allhic/plot.go
 * Path: /Users/bao/code/allhic
 * Created Date: Saturday, July 7th 2018, 1:33:37 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

// Plotter extracts a matrix of link counts and plot a heatmp
type Plotter struct {
	Anchor *Anchorer
}

// Run starts the plotting
func (r *Plotter) Run() {
	r.Anchor.Run()
}
