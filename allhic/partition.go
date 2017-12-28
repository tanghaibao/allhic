package allhic

import (
	"fmt"
	"os"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

// Partitioner is the bamfile into a matrix of link counts
type Partitioner struct {
	Bamfile string
}

// CountLinks provides the method to count the links
func (p Partitioner) CountLinks() {
	fmt.Println("In CountLinks()")
	fh, err := os.Open(p.Bamfile)
	if err != nil {
		panic(err)
	}

	br, err := bam.NewReader(fh, 0)
	if err != nil {
		panic(err)
	}
	defer br.Close()

	for _, rg := range br.Header().RGs() {
		fmt.Println(rg.Get(sam.Tag([2]byte{'S', 'M'})))
	}
}
