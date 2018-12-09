##############################################################################
###
### MST
###
### Reconstruct Random Main 5
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
load('RData/dataset.random.multiple.biopses.bulk.5.nodes.RData')

seed = 12345
available.cores = detectCores()

if (available.cores > 5) {
    cores = 5
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}


## Generate dataset for multiple.biopses 5

cat('result multiple.biopses.bulk 5\n')
result.random.multiple.biopses.bulk.5.nodes =
    expand.input(dataset.random.multiple.biopses.bulk.5.nodes,
                 seed,
                 cores,
                 pass.error.rates = FALSE)
save(result.random.multiple.biopses.bulk.5.nodes,
     file = "RData/result.random.multiple.biopses.bulk.5.nodes.RData")

### end of file -- reconstruct.random.main.5.R

