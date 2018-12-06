##############################################################################
###
### MST
###
### Plot Random Main 15
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
load('RData/dataset.random.multiple.biopses.15.nodes.RData')

seed = 12345
available.cores = detectCores()

if (available.cores > 5) {
    cores = 5
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}


## Generate dataset for multiple.biopses low

cat('result multiple.biopses 15\n')
result.random.multiple.biopses.15.nodes =
    expand.input(dataset.random.multiple.biopses.15.nodes,
                 seed,
                 cores,
                 pass.error.rates = FALSE)
save(result.random.multiple.biopses.15.nodes,
     file = "RData/result.random.multiple.biopses.15.nodes.RData")

### end of file -- reconstruct.random.main.15.R
