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
load('RData/dataset.random.single.cells.20.nodes.RData')

seed = 12345
available.cores = detectCores()

if (available.cores > 8) {
    cores = 8
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}

# generate dataset for single cells low
cat('result single cells 20\n')
result.random.single.cells.20.nodes = expand.input(dataset.random.single.cells.20.nodes, seed, cores)
save(result.random.single.cells.20.nodes, file="RData/result.random.single.cells.20.nodes.RData")
