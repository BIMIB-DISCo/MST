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
load('RData/dataset.random.multiple.biopses.bulk.15.nodes.RData')

seed = 12345
available.cores = detectCores()

if (available.cores > 5) {
    cores = 5
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}

# generate dataset for multiple.biopses.bulk low
cat('result multiple.biopses.bulk 15\n')
result.random.multiple.biopses.bulk.15.nodes = expand.input(dataset.random.multiple.biopses.bulk.15.nodes, seed, cores, pass.error.rates = FALSE)
save(result.random.multiple.biopses.bulk.15.nodes, file="RData/result.random.multiple.biopses.bulk.15.nodes.RData")
