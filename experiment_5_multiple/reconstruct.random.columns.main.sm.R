##############################################################################
###
### MST
###
### Reconstruct Random Columns Main SM
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

## Source the needed script

#  longjob -c "/afs/inf.ed.ac.uk/user/v/v1ldesa/download/R-3.2.5/bin/Rscript reconstruct.main.sh.R > high.log 2>&1"

library(parallel)

source('../reconstruct.run.R')
load('RData/dataset.multiple.biopses.random.columns.medium.RData')


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


## Generate dataset for multiple.biopses medium

cat('result multiple.biopses medium\n')
result.multiple.biopses.medium =
    expand.input(dataset.multiple.biopses.random.columns.medium,
                 seed,
                 cores,
                 pass.error.rates = FALSE)
save(result.multiple.biopses.medium,
     file = "RData/result.multiple.biopses.medium.RData")

### end of file -- reconstruct.random.columns.sm.R

