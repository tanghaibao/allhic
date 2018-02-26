/**
 * Filename: /Users/bao/code/allhic/agp.go
 * Path: /Users/bao/code/allhic
 * Created Date: Monday, February 26th 2018, 8:30:12 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

// AGPLine is a line in the AGP file
type AGPLine struct {
	object        string
	objectBeg     int
	objectEnd     int
	partNumber    int
	componentType byte
	isGap         bool
	strand        byte
	// As a gap
	gapLength       int
	gapType         string
	linkage         string
	linkageEvidence string
	// As a sequence chunk
	componentID  int
	componentBeg int
	componentEnd int
}
