/*
 *  extract_test.go
 *  allhic
 *
 *  Created by Haibao Tang on 12/11/19
 *  Copyright © 2019 Haibao Tang. All rights reserved.
 */

package allhic_test

import (
	"github.com/tanghaibao/allhic"
	"testing"
)

func TestCountSimplePattern(t *testing.T) {
	pattern := "GATC"
	seq := []byte("GATCGATCGATC")
	simplePattern := allhic.MakePattern(pattern)
	got := allhic.CountPattern(seq, simplePattern)
	expected := 3
	if got != expected {
		t.Errorf("CountPattern(#{seq}, #{pattern})=#{got}; want #{expected}")
	}
}

func TestCountRegexPattern(t *testing.T) {
	pattern := "GATCGATC,GANTGATC,GANTANTC,GATCANTC"
	seq := []byte("GATCGATCGGACTGATCGACCGATCACTCACGCTAAATGCAGAATCGATTATTC")
	regexPattern := allhic.MakePattern(pattern)
	got := allhic.CountPattern(seq, regexPattern)
	expected := 4
	if got != expected {
		t.Errorf("CountPattern(#{seq}, #{pattern})=#{got}; want #{expected}")
	}
}
