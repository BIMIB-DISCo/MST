##############################################################################
###
### MST
###
### Reconstruct Main SCs Mini 05
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
load('RData/mini_dataset_05.RData')


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


## Generate dataset for single cells high

cat('result single cells mini\n')
result.single.cells.mini.05 = expand.input(mini_dataset_05, seed, cores)
save(result.single.cells.mini.05,
     file = "RData/result.single.cells.mini.05.RData")

### end of file -- reconstruct.main.smini.05.R
