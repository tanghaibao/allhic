import future
import parseutils
import strutils
import tables


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


proc initCLMFile*(name: string, clmfile: string): CLMFile =
  new result
  result.name = name
  result.clmfile = clmfile
  result.idsfile = clmfile.rsplit('.', maxsplit=1)[0] & ".ids"
  result.tig_to_size = initOrderedTable[string, int]()
  result.contacts = initTable[(string, string), int]()


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
