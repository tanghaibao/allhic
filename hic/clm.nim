import future
import math
import os
import parseutils
import sets
import strutils
import tables
import logging
import "matrix"

var logger = newConsoleLogger()
logger.addHandler()


##
## CLM file (modified) has the following format:
##    tig00046211+ tig00063795+       1       53173
##    tig00046211+ tig00063795-       1       116050
##    tig00046211- tig00063795+       1       71155
##    tig00046211- tig00063795-       1       134032
##    tig00030676+ tig00077819+       7       136407 87625 87625 106905 102218 169660 169660
##    tig00030676+ tig00077819-       7       126178 152952 152952 35680 118923 98367 98367
##    tig00030676- tig00077819+       7       118651 91877 91877 209149 125906 146462 146462
##    tig00030676- tig00077819-       7       108422 157204 157204 137924 142611 75169 75169
##
type
  CLMFile* = ref object
    name*: string
    clmfile*: string
    idsfile*: string
    tig_to_size*: OrderedTable[string, int]
    contacts*: Table[(string, string), int]
    active*: OrderedSet[string]


proc initCLMFile*(name: string, clmfile: string): CLMFile =
  new result
  result.name = name
  result.clmfile = clmfile
  result.idsfile = clmfile.rsplit('.', maxsplit=1)[0] & ".ids"
  result.tig_to_size = initOrderedTable[string, int]()
  result.contacts = initTable[(string, string), int]()
  result.active = initOrderedSet[string]()


##
## IDS file has a list of contigs that need to be ordered. 'recover',
## keyword, if available in the third column, is less confident.
##        tig00015093     46912
##        tig00035238     46779   recover
##        tig00030900     119291
##
proc parse_ids*(this: CLMFile, skiprecover: bool) =
  for line in lines this.idsfile:
    let
      atoms = line.split()
      tig = atoms[0]
      size = parseInt(atoms[1])
    this.tig_to_size[tig] = size
    this.active.incl(tig)


proc parse_clm*(this: CLMFile) =
  for line in lines this.clmfile:
    let atoms = line.strip().split('\t')
    doAssert len(atoms) == 3, "Malformed line `$#`".format(atoms)
    let
      abtig = atoms[0].split()
      links = parseInt(atoms[1])
      dists = lc[parseInt(x) | (x <- atoms[2].split()), int]
      atig = abtig[0]
      btig = abtig[1]
      at = atig[0..^2]
      ao = atig[^1]
      bt = btig[0..^2]
      bo = btig[^1]

    if at notin this.tig_to_size:
      continue
    if bt notin this.tig_to_size:
      continue

    this.contacts[(at, bt)] = links


proc N*(this: CLMFile): int =
  len(this.active)


proc tig_to_idx*(this: CLMFile): Table[string, int] =
  result = initTable[string, int]()
  var i = 0
  for key in this.active:
    result[key] = i
    inc i


##
## Contact frequency matrix. Each cell contains how many inter-contig links
## between i-th and j-th contigs.
##
proc M*(this: CLMFile): Matrix[int] =
  let
    N = this.N
    tig_to_idx = this.tig_to_idx

  result = newMatrix[int](N, N)
  for abt, links in this.contacts.pairs():
    let (at, bt) = abt
    if at notin tig_to_idx:
      continue
    if bt notin tig_to_idx:
      continue
    let
      ai = tig_to_idx[at]
      bi = tig_to_idx[bt]
    result[ai, bi] = links
    result[bi, ai] = links


proc cumSum*[T](x: openArray[T]): seq[T] =
  ##
  ## Cumulative sum for each element of ``x``
  ##
  ## ``cumSum(@[1,2,3,4])`` produces ``@[1,3,6,10]``
  ##
  result = newSeq[T](x.len)
  var cp = T(0)
  for i in 0..<x.len:
    cp = cp + x[i]
    result[i] = cp


proc score_evaluate*(tour: seq[int], tour_sizes: seq[int], tour_M: Matrix[int]): int =
  const
    LIMIT = 10_000_000

  var sizes_cum = newSeq[int](tour.len)
  var cp = 0
  var mid = 0
  for i, t in tour.pairs():
    sizes_cum[i] = cp + tour_sizes[t] div 2
    cp = cp + tour_sizes[t]

  var
    s = 0.0
    a, b, ia, ib, dist, links: int

  let size = tour.len
  for ia in 0..<size:
    a = tour[ia]
    for ib in (ia + 1)..<size:
      b = tour[ib]
      links = tour_M[a, b]
      if links == 0:
        continue
      dist = sizes_cum[ib] - sizes_cum[ia]
      if dist > LIMIT:
        break
      s += links / dist

  echo s
  echo tour
  echo tour_sizes
  echo sizes_cum


proc active_sizes*(this: CLMFile): seq[int] =
  lc[this.tig_to_size[x] | (x <- this.active), int]


proc report_active*(this: CLMFile) =
  debug("Active contigs: $# (length=$#)".format(this.N, this.active_sizes.sum()))


##
## Select contigs in the current partition. This is the setup phase of the
## algorithm, and supports two modes:
##  - "de novo": This is useful at the start of a new run where no tours are
##  available. We select the strong contigs that have significant number
##  of links to other contigs in the partition. We build a histogram of
##  link density (# links per bp) and remove the contigs that appear to be
##  outliers. The orientations are derived from the matrix decomposition
##  of the pairwise strandedness matrix O.
##  - "hotstart": This is useful when there was a past run, with a given
##  tourfile. In this case, the active contig list and orientations are
##  derived from the last tour in the file.
##
proc activate*(this: CLMFile, tourfile = "", minsize = 10000) =
  var tf = tourfile
  if existsFile(tourfile):
    tf = ""

  if true:
    this.report_active()

  debug("Remove contigs with size < $#".format(minsize))
  this.active = lc[x | (x <- this.active, this.tig_to_size[x] >= minsize),
                string].toOrderedSet()
  this.report_active()


proc parse*(this: CLMFile) =
  this.parse_ids(true)
  this.parse_clm()
  this.activate()

  let M = this.M
  let N = this.N
  let tour = lc[x | (x <- 0..N - 1), int]
  let tour_sizes = this.active_sizes()
  discard score_evaluate(tour, tour_sizes, M)
