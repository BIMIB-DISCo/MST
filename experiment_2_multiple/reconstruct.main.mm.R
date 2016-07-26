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

#  longjob -c "/afs/inf.ed.ac.uk/user/v/v1ldesa/download/R-3.2.5/bin/Rscript reconstruct.main.sh.R > high.log 2>&1"

library(parallel)

source('../reconstruct.run.R')
load('RData/dataset.multiple.biopses.medium.RData')

# setting of the experiments
seed = 12345

available.cores = detectCores()

if (available.cores > 5) {
    cores = 5
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}

# generate dataset for multiple.biopses medium
cat('result multiple.biopses medium\n')
result.multiple.biopses.medium = expand.input(dataset.multiple.biopses.medium, seed, cores, pass.error.rates = FALSE)
save(result.multiple.biopses.medium, file="RData/result.multiple.biopses.medium.RData")
