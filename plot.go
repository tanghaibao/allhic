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

// Plotter extracts a matrix of link counts and plot a heatmp
type Plotter struct {
	Anchor *Anchorer
}

// Run starts the plotting
func (r *Plotter) Run() error {
	m := r.Anchor
	err := m.ExtractInterContigLinks()
	if err != nil {
		return err
	}
	m.parseTourFile(m.Tourfile)
	err = m.printTour(os.Stdout, "PLOTTER")
	if err != nil {
		return err
	}
	// Serialize to disk for plotting
	m.makeContigStarts()
	err = m.serialize(250000, "genome.json", "data.npy")
	if err != nil {
		return err
	}
	err = r.host()
	if err != nil {
		return err
	}
	log.Notice("Success")
	return nil
}

// Host plot.html
func (r *Plotter) host() error {
	box := packr.NewBox("./templates")
	port := 3000
	f, _ := os.Create("index.html")
	s, _ := box.FindString("index.html")
	_, err := f.WriteString(s)
	if err != nil {
		return err
	}
	err = f.Sync()
	if err != nil {
		return err
	}

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
	return f.Close()
}
