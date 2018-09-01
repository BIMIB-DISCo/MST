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

# this is based on experiment 3 single cell

# source the needed script
library(TRONCO)
library(parallel)

library(igraph)

if (!dir.exists('RData')) {
    dir.create('RData')
}


source('../generate.dataset.R')
source('../generate.scite.input.R')

# setting of the experiments
seed = 12345
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
scite.sd = 0

# setting for single cell
sample_sizes_single_cells = c(100, 500, 1000, 5000, 10000)
e_pos_single_cells = c(0.005)
e_neg_single_cells = c(0.050)

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
clusterExport(cl, c('sample_sizes_single_cells',
    'e_pos_single_cells',
    'e_neg_single_cells'))
clusterSetRNGStream(cl, iseed = seed)

cat('Using', cores, 'cores via "parallel" \n')

# generate dataset for random single cells of 20 nodes
cat('dataset random single cells of 20 nodes\n')
dataset.random.single.cells.20.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells(type="random",
        samples_num = sample_sizes_single_cells,
        e_pos = e_pos_single_cells,
        e_neg = e_neg_single_cells,
        nodes = 20,
        min_significance = 0.60,
        max_significance = 0.90,
        samples_significance = 0.001)
})
save(dataset.random.single.cells.20.nodes, file="RData/dataset.random.single.cells.20.nodes.RData")

stopCluster(cl)

#print('scite single')
#create.scite.input(dataset.random.single.cells.20.nodes, 'single', 'random_20', scite.sd)
