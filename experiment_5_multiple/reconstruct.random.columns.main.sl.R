##############################################################################
###
### MST
###
### Reconstruct Random Columns Main SL
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
load('RData/dataset.multiple.biopses.random.columns.low.RData')


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


## Generate dataset for multiple.biopses low

cat('result multiple.biopses low\n')
result.multiple.biopses.low =
    expand.input(dataset.multiple.biopses.random.columns.low,
                 seed,
                 cores,
                 pass.error.rates = FALSE)
save(result.multiple.biopses.low,
     file = "RData/result.multiple.biopses.low.RData")

### end of file -- reconstruct.random.columns.main.sl.R
