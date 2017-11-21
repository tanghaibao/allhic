# Package

version       = "0.7.11"
author        = "Haibao Tang"
description   = "Genome scaffolding based on HiC data"
license       = "BSD"

# Dependencies

requires "nim >= 0.17.2"

task run, "build":
    exec "nim c -d:release -r allhic.nim"
