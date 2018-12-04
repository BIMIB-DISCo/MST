##############################################################################
###
### MST
###
### Reconstruct Main Convergent
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

## Source the needed script

library(parallel)

source('../reconstruct.run.R')
load('RData/dataset.multiple.biopses.convergent.RData')

## Setting of the experiments

seed = 12345

available.cores = detectCores()

if (available.cores > 8) {
    cores = 8
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}


## Generate dataset for single cells low

cat('result single cells low\n')
result.multiple.biopses.convergent = expand.input(dataset.multiple.biopses.convergent, seed, cores)
save(result.multiple.biopses.convergent, file = "RData/result.multiple.biopses.convergent.RData")

### end of file -- reconstruct.main.convergent
