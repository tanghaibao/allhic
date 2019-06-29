/*
 *  clm_test.go
 *  allhic
 *
 *  Created by Haibao Tang on 06/28/19
 *  Copyright © 2019 Haibao Tang. All rights reserved.
 */

package allhic_test

import (
	"path"
	"testing"

	"github.com/tanghaibao/allhic"
)

func TestParseRECountsFile(t *testing.T) {
	reCountsFile := allhic.RECountsFile{
		Filename: path.Join("tests", "test.counts_RE.txt"),
	}
	reCountsFile.ParseRecords()
	expectedNumRecords := 2
	if len(reCountsFile.Records) != expectedNumRecords {
		t.Fatalf("Expected %d records, got %d records", expectedNumRecords, len(reCountsFile.Records))
	}
}
