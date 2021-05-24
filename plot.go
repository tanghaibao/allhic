/*
 * Filename: /Users/bao/code/allhic/plot.go
 * Path: /Users/bao/code/allhic
 * Created Date: Saturday, July 7th 2018, 1:33:37 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"net/http"
	"os"
	"strconv"

	"github.com/gobuffalo/packr"
)

// Plotter extracts a matrix of link counts and plot a heatmap
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
	r.host()

	// printPaths(paths)
	log.Notice("Success")
}

// Host plot.html
func (r *Plotter) host() {
	box := packr.NewBox("./templates")
	port := 3000
	f, _ := os.Create("index.html")
	s, _ := box.FindString("index.html")
	_, _ = f.WriteString(s)
	_ = f.Sync()

	http.Handle("/", http.FileServer(http.Dir(".")))

	for {
		log.Noticef("Serving on localhost:%d ...", port)
		if err := http.ListenAndServe(":"+strconv.Itoa(port), nil); err != nil {
			log.Debug(err)
			port++
		} else {
			break
		}
	}
	_ = f.Close()
}
