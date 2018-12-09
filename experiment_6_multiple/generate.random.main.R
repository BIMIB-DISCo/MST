##############################################################################
###
### MST
###
### Generate Random Main
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

## Source the needed script

library(TRONCO)
library(parallel)

library(igraph)

if (!dir.exists('RData')) {
    dir.create('RData')
}


source('../generate.dataset.R')
source('../generate.scite.input.R')
source('../generate.sample.polyclonal.trees.R')


## Setting of the experiments

seed = 12345
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
scite.sd = 0

## Setting for single cell

wild_type_rate = 0.1
sample_sizes_multiple_biopses = c(5, 7, 10, 20, 50)
e_pos_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200)
e_neg_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200)

available.cores = detectCores()

if (available.cores > 8) {
    cores = 8
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}

cl = makeCluster(cores)
clusterEvalQ(cl, library(igraph))
clusterEvalQ(cl, source('../generate.dataset.R'))
clusterEvalQ(cl, source('../generate.sample.polyclonal.trees.R'))
clusterExport(cl, c('sample_sizes_multiple_biopses',
    'e_pos_multiple_biopses',
    'e_neg_multiple_biopses',
    'wild_type_rate'))
clusterSetRNGStream(cl, iseed = seed)

cat('Using', cores, 'cores via "parallel" \n')

## Generate dataset for random single cells of 5 nodes

cat('dataset random single cells forest of 20 nodes\n')
dataset.random.multiple.biopses.forest.nodes =
    parSapply(cl,
              my_experiments,
              function(x) {
                  generate.dataset.multiple.biopses(type = "random_forest",
                                                    samples_num = sample_sizes_multiple_biopses,
                                                    e_pos = e_pos_multiple_biopses,
                                                    e_neg = e_neg_multiple_biopses,
                                                    nodes = 20,
                                                    min_significance = 0.70,
                                                    max_significance = 0.90,
                                                    samples_significance = 0.001,
                                                    wild_type = wild_type_rate)
              })
save(dataset.random.multiple.biopses.forest.nodes,
     file = "RData/dataset.random.multiple.biopses.forest.nodes.RData")


stopCluster(cl)

print('scite multiple')
create.scite.input(dataset.random.multiple.biopses.forest.nodes,
                   'multiple',
                   'random_forest',
                   scite.sd)

### end of file -- generate.random.main.R
