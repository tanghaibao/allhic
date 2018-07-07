/*
 * Filename: /Users/bao/code/allhic/model.go
 * Path: /Users/bao/code/allhic
 * Created Date: Friday, July 6th 2018, 10:47:29 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"bufio"
	"fmt"
	"math"
	"os"
)

// LinkDensityModel is a power-law model Y = A * X ^ B, stores co-efficients
// this density than needs to multiply C - X to make it a probability distribution
// where C is chromosome length
type LinkDensityModel struct {
	A, B        float64
	binStarts   []int
	binNorms    []int
	nLinks      []int
	linkDensity []float64
}

// ********* Calculation of link distribution model ************

// NewLinkDensityModel makes an empty link distribution ready to be filled in
func NewLinkDensityModel() *LinkDensityModel {
	return &LinkDensityModel{
		binStarts:   []int{},
		binNorms:    make([]int, nBins),
		nLinks:      make([]int, nBins),
		linkDensity: make([]float64, nBins),
	}
}

// writeDistribution writes the link size distribution to file
func (r *LinkDensityModel) writeDistribution(outfile string) {
	f, _ := os.Create(outfile)
	w := bufio.NewWriter(f)
	defer f.Close()

	fmt.Fprintf(w, "#Bin\tBinStart\tBinSize\tNumLinks\tTotalSize\tLinkDensity\n")
	for i := 0; i < nBins; i++ {
		fmt.Fprintf(w, "%d\t%d\t%d\t%d\t%d\t%.4g\n",
			i, r.binStarts[i], r.BinSize(i), r.nLinks[i], r.binNorms[i], r.linkDensity[i])
	}

	w.Flush()
	log.Noticef("Link size distribution written to `%s`", outfile)
}

// linkBin takes a link distance and convert to a binID
func (r *LinkDensityModel) linkBin(dist int) int {
	if dist < MinLinkDist {
		return -1
	}
	distOverMin := dist / MinLinkDist
	log2i := uintLog2(uint(distOverMin))
	log2f := uintLog2Frac(float64(dist) / float64(int(MinLinkDist)<<log2i))
	return int(16*log2i + log2f)
}

// BinSize returns the size of each bin
func (r *LinkDensityModel) BinSize(i int) int {
	return r.binStarts[i+1] - r.binStarts[i]
}

// makeNorms computes the normalization size for each bin
func (r *LinkDensityModel) makeNorms(contigSizes []int) {
	for _, size := range contigSizes {
		for j := 0; j < nBins; j++ {
			z := size - r.binStarts[j]
			if z < 0 {
				break
			}
			r.binNorms[j] += z
		}
	}
}

// makeBins makes geometric bins and count links that fall in each bin
// This heavily borrows the method form LACHESIS
// https://github.com/shendurelab/LACHESIS/blob/master/src/LinkSizeExtracter.cc
func (r *LinkDensityModel) makeBins() {
	// Step 1: make geometrically sized bins
	for i := 0; 16*i <= nBins; i++ {
		jpower := 1.0
		for j := 0; j < 16 && 16*i+j <= nBins; j++ {
			binStart := MinLinkDist << uint(i)
			r.binStarts = append(r.binStarts, int(float64(binStart)*jpower))
			jpower *= GeometricBinSize
		}
	}
}

// countBinDensities counts the links in each bin and divide by the biniesorm
func (r *LinkDensityModel) countBinDensities(contigs []*ContigInfo) {
	maxLinkDist := math.MinInt32
	for _, contig := range contigs {
		for _, link := range contig.links {
			if link > maxLinkDist {
				maxLinkDist = link
			}
		}
	}
	// Step 2: calculate assayable sequence length
	// Find the length of assayable intra-contig sequence in each bin
	intraContigLinkRange := math.Log2(float64(maxLinkDist) / float64(MinLinkDist))
	nIntraContigBins := int(math.Ceil(intraContigLinkRange * 16))

	// Step 3: loop through all links and tabulate the counts
	for _, contig := range contigs {
		for _, link := range contig.links {
			bin := r.linkBin(link)
			if bin == -1 {
				continue
			}
			r.nLinks[bin]++
		}
	}

	// Step 4: normalize to calculate link density
	for i := 0; i < nIntraContigBins; i++ {
		r.linkDensity[i] = float64(r.nLinks[i]) / float64(r.binNorms[i]) / float64(r.BinSize(i))
	}

	// Step 5: in LACHESIS, we assume the distribution approximates 1/x
	// for large x. This is not accurate. We should instead infer a power law
	// distribution.
	topBin := nIntraContigBins - 1
	nTopLinks := 0
	nTopLinksNeeded := len(contigs[0].links) / 100
	for ; nTopLinks < nTopLinksNeeded; topBin-- {
		nTopLinks += r.nLinks[topBin]
	}

	Xs := []int{}
	Ys := []float64{}
	for i := 0; i < topBin; i++ {
		Xs = append(Xs, r.binStarts[i])
		Ys = append(Ys, r.linkDensity[i])
	}
	r.fitPowerLaw(Xs, Ys)

	// Overwrite the values of last few bins, or a bin with na values
	for i := 0; i < nBins; i++ {
		if r.linkDensity[i] == 0 || i >= topBin {
			r.linkDensity[i] = r.transformPowerLaw(r.binStarts[i])
		}
	}
}

// fitPowerLaw fits power law distribution
// See reference: http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
// Assumes the form Y = A * X^B, returns (A, B), the coefficients
func (r *LinkDensityModel) fitPowerLaw(Xs []int, Ys []float64) {
	SumLogXLogY, SumLogXLogX, SumLogX, SumLogY := 0.0, 0.0, 0.0, 0.0
	n := len(Xs)
	for i := 0; i < n; i++ {
		logXs, logYs := math.Log(float64(Xs[i])), math.Log(Ys[i])
		SumLogXLogY += logXs * logYs
		SumLogXLogX += logXs * logXs
		SumLogX += logXs
		SumLogY += logYs
	}

	B := (float64(n)*SumLogXLogY - SumLogX*SumLogY) / (float64(n)*SumLogXLogX - SumLogX*SumLogX)
	A := math.Exp((SumLogY - B*SumLogX) / float64(n))
	r.A, r.B = A, B

	log.Noticef("Power law Y = %.4f * X ^ %.4f", A, B)
}

// transformPowerLaw interpolate probability value given a link size
func (r *LinkDensityModel) transformPowerLaw(X int) float64 {
	return r.A * math.Pow(float64(X), r.B)
}

// transformLogProb calculates the probability given a link size
func (r *LinkDensityModel) tranformLogProb(X int) float64 {
	// The following two version have subtle differences, first one is more accurate, but
	// in reality the difference seems to be neglible
	// return math.Log(float64(r.Seqsize-X)) + r.model.B*math.Log(float64(X))
	return r.B * math.Log(float64(X))
}
