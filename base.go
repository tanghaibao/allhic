/**
 * Filename: /Users/bao/code/allhic/allhic/base.go
 * Path: /Users/bao/code/allhic/allhic
 * Created Date: Tuesday, January 2nd 2018, 8:07:22 pm
 * Author: bao
 *
 * Copyright (c) 2018 Haibao Tang
 */

package allhic

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"math"
	"os"
	"path"
	"sort"
	"strings"

	"github.com/op/go-logging"
)

const (
	// Version is the current version of ALLHIC
	Version = "0.9.13"
	// LB is lower bound for GoldenArray
	LB = 18
	// UB is upper bound for GoldenArray
	UB = 29
	// BB is span for GoldenArray
	BB = UB - LB + 1
	// PHI is natural log of golden ratio
	PHI = 0.4812118250596684 // math.Log(1.61803398875)
	// OUTLIERTHRESHOLD is how many deviation from MAD
	OUTLIERTHRESHOLD = 3.5
	// MINSIZE is the minimum size cutoff for tig to be considered
	MINSIZE = 10000
	// GeometricBinSize is the max/min ratio for each bin
	GeometricBinSize = 1.0442737824274138403219664787399
	// MinLinkDist is the minimum link distance we care about
	MinLinkDist = 1 << 11

	/* extract */

	// DefaultRE is the default restriction site used
	DefaultRE = "GATC"
	// MinLinks is the minimum number of links between contig pair to consider
	MinLinks = 3

	// MaxLinkDist is the maximum link distance we care about
	MaxLinkDist = 1 << 27
	// BigNorm is a big integer multiplier so we don't have to mess with float64
	BigNorm = int64(1000000000000)
	// MinAvgLinkage is the minimum cutoff for merging clusters
	MinAvgLinkage = 0
	// LinkDist specifies to maximum size of the links going over a certain position
	LinkDist = int64(1000000)

	/* optimize */

	// Seed is the random seed
	Seed = 42
	// Npop is the population size used in GA
	Npop = 100
	// Ngen is the number of generations for convergence
	Ngen = 5000
	// MutaProb is the mutation probability in GA
	MutaProb = 0.2

	// *** The following parameters are modeled after LACHESIS ***

	// MinREs is the minimum number of RE sites in a contig to be clustered (CLUSTER_MIN_RE_SITES)
	MinREs = 10
	// MaxLinkDensity is the density threshold before marking a contig as 'repetitive' (CLUSTER_MAX_LINK_DENSITY)
	MaxLinkDensity = 2
	// NonInformativeRatio is the cutoff for recovering skipped contigs back into the clusters (CLUSTER_NON-INFORMATIVE_RATIO)
	NonInformativeRatio = 3

	// *** CSV headers ***

	// REHeader is the first line in the RE counts file
	REHeader = "#Contig\tRECounts\tLength\n"

	// PairsFileHeader is the first line in the pairs.txt file
	PairsFileHeader = "#X\tY\tContig1\tContig2\tRE1\tRE2\tObservedLinks\tExpectedLinksIfAdjacent\tLabel\n"

	// DistributionHeader is the first line in the distribution.txt file
	DistributionHeader = "#Bin\tBinStart\tBinSize\tNumLinks\tTotalSize\tLinkDensity\n"

	// PostProbHeader is the first line in the postprob file
	PostProbHeader = "#SeqID\tStart\tEnd\tContig\tPostProb\n"
)

// GArray contains golden array of size BB
type GArray [BB]int

// GR is a precomputed list of exponents of golden ratio phi
var GR = [...]int{5778, 9349, 15127, 24476,
	39603, 64079, 103682, 167761,
	271443, 439204, 710647, 1149851}

var log = logging.MustGetLogger("allhic")
var format = logging.MustStringFormatter(
	`%{color}%{time:15:04:05} %{shortfunc} | %{level:.6s} %{color:reset} %{message}`,
)

// Backend is the default stderr output
var Backend = logging.NewLogBackend(os.Stderr, "", 0)

// BackendFormatter contains the fancy debug formatter
var BackendFormatter = logging.NewBackendFormatter(Backend, format)

// ErrorAbort logs an error message and then exit with retcode of 1
func ErrorAbort(err error) {
	if err != nil {
		log.Errorf("%s", err)
		os.Exit(1)
	}
}

// RemoveExt returns the substring minus the extension
func RemoveExt(filename string) string {
	return strings.TrimSuffix(filename, path.Ext(filename))
}

// Round makes a round number
func Round(input float64) float64 {
	if input < 0 {
		return math.Ceil(input - 0.5)
	}
	return math.Floor(input + 0.5)
}

// SumLog returns the kernel of sum of log likelihood
func SumLog(a []int) float64 {
	sum := 0.0
	for _, val := range a {
		sum += math.Log(float64(val))
	}
	return sum
}

// GoldenArray is given list of ints, we aggregate similar values so that it becomes an
// array of multiples of phi, where phi is the golden ratio.
//
// phi ^ 18 = 5778
// phi ^ 29 = 1149851
//
// So the array of counts go between 843 to 788196. One triva is that the
// exponents of phi gets closer to integers as N grows. See interesting
// discussion here:
// <https://www.johndcook.com/blog/2017/03/22/golden-powers-are-nearly-integers/>
func GoldenArray(a []int) (counts GArray) {
	for _, x := range a {
		c := int(Round(math.Log(float64(x)) / PHI))
		if c < LB {
			c = LB
		} else if c > UB {
			c = UB
		}
		counts[c-LB]++
	}
	return
}

// abs gets the absolute value of an int
func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// min gets the minimum for two ints
func min(x, y int) int {
	if x < y {
		return x
	}
	return y
}

// minInt64 gets the minimum for two int64
func minInt64(x, y int64) int64 {
	if x < y {
		return x
	}
	return y
}

// max gets the maximum for two ints
func max(x, y int) int {
	if x > y {
		return x
	}
	return y
}

// sumf gets the sum for an int slice
func sumf(a []float64) float64 {
	ans := 0.0
	for _, x := range a {
		ans += x
	}
	return ans
}

// unique returns a distinct slice of ints
func unique(a []int) []int {
	keys := make(map[int]bool)
	list := make([]int, 0)
	for _, entry := range a {
		if _, value := keys[entry]; !value {
			keys[entry] = true
			list = append(list, entry)
		}
	}
	sort.Ints(list)
	return list
}

// arrayToString print delim-separated int slice
func arrayToString(a []int, delim string) string {
	sort.Ints(a)
	return strings.Trim(strings.Replace(fmt.Sprint(a), " ", delim, -1), "[]")
}

// median gets the median value of an array
func median(s []float64) float64 {
	// Make a sorted copy
	numbers := make([]float64, len(s))
	copy(numbers, s)
	sort.Float64s(numbers)

	middle := len(numbers) / 2
	result := numbers[middle]
	if len(numbers)%2 == 0 {
		result = (result + numbers[middle-1]) / 2
	}
	return result
}

// L50 returns the sequence length L where half of the genome is covered in contigs
// of length >= L50
func L50(lengths []int64) int64 {
	sortInt64s(lengths)
	total := int64(0)
	for _, length := range lengths {
		total += length
	}
	halfTotal := total / 2
	cumsize := int64(0)
	i := len(lengths) - 1
	for ; i >= 0; i-- {
		cumsize += lengths[i]
		if cumsize > halfTotal {
			break
		}
	}
	return lengths[i]
}

// OutlierCutoff implements Iglewicz and Hoaglin's robust, returns the cutoff values -
// lower bound and upper bound.
func OutlierCutoff(a []float64) (float64, float64) {
	M := median(a)
	D := make([]float64, len(a))
	for i := 0; i < len(a); i++ {
		D[i] = math.Abs(a[i] - M)
	}
	MAD := median(D)
	C := OUTLIERTHRESHOLD / .67449 * MAD
	return M - C, M + C
}

// Make2DSlice allocates a 2D matrix with shape (m, n)
func Make2DSlice(m, n int) [][]int {
	P := make([][]int, m)
	for i := 0; i < m; i++ {
		P[i] = make([]int, n)
	}
	return P
}

// Make2DSliceInt64 allocates a 2D int64 matrix with shape (m, n)
func Make2DSliceInt64(m, n int) [][]int64 {
	P := make([][]int64, m)
	for i := 0; i < m; i++ {
		P[i] = make([]int64, n)
	}
	return P
}

// Make2DGArraySlice allocates a 2D matrix with shape (m, n)
func Make2DGArraySlice(m, n int) [][]GArray {
	P := make([][]GArray, m)
	for i := 0; i < m; i++ {
		P[i] = make([]GArray, n)
	}
	return P
}

// Percentage prints a human readable message of the percentage
func Percentage(a, b int) string {
	return fmt.Sprintf("%d of %d (%.1f %%)", a, b, float64(a)*100./float64(b))
}

// ReadCSVLines parses all the csv lines into 2D array of tokens
func ReadCSVLines(filename string) [][]string {
	log.Noticef("Parse csvfile `%s`", filename)
	fh := mustOpen(filename)

	var data [][]string
	r := csv.NewReader(bufio.NewReader(fh))
	r.Comma = '\t'
	for i := 0; ; i++ {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatal(err)
		}
		if i == 0 {
			continue // Skip header
		}
		data = append(data, rec)
	}

	_ = fh.Close()
	return data
}

// sortInt64s sorts a slice of int64
func sortInt64s(a []int64) {
	sort.Slice(a, func(i, j int) bool {
		return a[i] < a[j]
	})
}

// searchInt64 searches the position within a sorted int64 slice
func searchInt64(a []int64, x int64) int {
	return sort.Search(len(a), func(i int) bool { return a[i] >= x })
}

// mustExist panics when a file is not found
func mustExist(filename string) {
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		log.Fatal(err)
	}
}

// mustOpen wraps os.Open but panics if file not found
func mustOpen(filename string) *os.File {
	f, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	return f
}
