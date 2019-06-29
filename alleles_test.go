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

func TestParsePafFile(t *testing.T) {
	path := filepath.Join("tests", "test.paf")
	alleler := allhic.Alleler{PafFile: path}
	paf := alleler.Run()
	if len(paf.Records) != 10 {
		t.Fatal("Number of records != 10")
	}
	cmValue, ok := paf.Records[0].Tags["cm"].(int)
	if !ok {
		t.Fatal("Cannot find the cm tag in the first PAF record")
	}
	expectedCmValue := 5277
	if cmValue != expectedCmValue {
		t.Fatalf("The first record is expected to have cm = %d, got %d", expectedCmValue, cmValue)
	}
}
