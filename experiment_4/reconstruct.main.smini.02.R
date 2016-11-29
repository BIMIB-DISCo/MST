##################################################################################
#                                                                                #
# MST                                                                            #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Giulio Caravagna, Luca De Sano, Daniele Ramazzotti         #
# email: tronco@disco.unimib.it                                                  #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################

# source the needed script

library(parallel)

source('../reconstruct.run.R')
load('RData/mini_dataset_02.RData')

# setting of the experiments
seed = 12345

available.cores = detectCores()

if (available.cores > 8) {
    cores = 8
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}

# generate dataset for single cells high
cat('result single cells mini\n')
result.single.cells.mini.02 = expand.input(mini_dataset_02, seed, cores)
save(result.single.cells.mini.02, file="RData/result.single.cells.mini.02.RData")
