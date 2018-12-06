##############################################################################
###
### MST
###
### Reconstruct Main MH
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
load('RData/dataset.multiple.biopses.high.RData')


## Setting of the experiments

seed = 12345

available.cores = detectCores()

if (available.cores > 5) {
    cores = 5
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}


## Generate dataset for multiple.biopses high

cat('result multiple.biopses high\n')
result.multiple.biopses.high =
    expand.input(dataset.multiple.biopses.high,
                 seed,
                 cores,
                 pass.error.rates = FALSE)
save(result.multiple.biopses.high,
     file = "RData/result.multiple.biopses.high.RData")


### end of file -- reconstruct.main.mh.R
