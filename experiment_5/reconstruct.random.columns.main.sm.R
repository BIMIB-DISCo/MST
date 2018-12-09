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
load('RData/dataset.single.cells.random.columns.medium.RData')


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


## Generate dataset for single cells medium

cat('result single cells medium\n')
result.single.cells.medium =
    expand.input(dataset.single.cells.random.columns.medium, seed, cores)
save(result.single.cells.medium,
     file = "RData/result.single.cells.medium.RData")

### end of file -- reconstruct.random.columns.main.sm.R
