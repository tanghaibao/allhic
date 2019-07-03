/*
 *  alleles_test.go
 *  allhic
 *
 *  Created by Haibao Tang on 06/28/19
 *  Copyright © 2019 Haibao Tang. All rights reserved.
 */

package allhic_test

import (
	"path/filepath"
	"testing"

	"github.com/tanghaibao/allhic"
)

// setupAlleler reads in test.paf and return the object for testing
func setupAlleler() allhic.Alleler {
	pafFile := filepath.Join("tests", "test.paf")
	reFile := filepath.Join("tests", "test.counts_RE.txt")
	alleler := allhic.Alleler{PafFile: pafFile, ReFile: reFile}
	alleler.Run()
	return alleler
}

func TestParsePafFile(t *testing.T) {
	alleler := setupAlleler()
	expectedNumRecords := 10
	if len(alleler.Paf.Records) != expectedNumRecords {
		t.Fatalf("Expected %d records, got %d", expectedNumRecords, len(alleler.Paf.Records))
	}
	cmValue, ok := alleler.Paf.Records[0].Tags["cm"].(int)
	if !ok {
		t.Fatal("Cannot find the cm tag in the first PAF record")
	}
	expectedCmValue := 5277
	if cmValue != expectedCmValue {
		t.Fatalf("The first record is expected to have cm %d, got %d", expectedCmValue, cmValue)
	}
	expectedLength := 135917
	if gotLength := alleler.ReCounts.Records[0].Length; gotLength != expectedLength {
		t.Fatalf("The first record is expected to have length %d, got %d", expectedLength, gotLength)
	}
}
