# Package

version       = "0.7.11"
author        = "Haibao Tang"
description   = "Genome scaffolding based on HiC data"
license       = "BSD"

# Dependencies

requires "nim >= 0.17.2"

task run, "run":
    exec "nim c -d:release allhic.nim " &
     "&& LD_LIBRARY_PATH=./htslib ./allhic partition"
