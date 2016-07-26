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

# setting for multiple biopses
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

# generate dataset for random multiple.biopses of 5 nodes
cat('dataset random multiple.biopses of 5 nodes\n')
dataset.random.multiple.biopses.5.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses(type="random",
        samples_num = sample_sizes_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate,
        nodes = 5,
        min_significance = 0.60,
        max_significance = 0.90,
        samples_significance = 0.001)
})
save(dataset.random.multiple.biopses.5.nodes, file="RData/dataset.random.multiple.biopses.5.nodes.RData")

# generate dataset for random multiple.biopses of 10 nodes
cat('dataset random multiple.biopses of 10 nodes\n')
dataset.random.multiple.biopses.10.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses(type="random",
        samples_num = sample_sizes_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate,
        nodes = 10,
        min_significance = 0.60,
        max_significance = 0.90,
        samples_significance = 0.001)
})
save(dataset.random.multiple.biopses.10.nodes, file="RData/dataset.random.multiple.biopses.10.nodes.RData")

# generate dataset for random multiple.biopses of 15 nodes
cat('dataset random multiple.biopses of 15 nodes\n')
dataset.random.multiple.biopses.15.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses(type="random",
        samples_num = sample_sizes_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate,
        nodes = 15,
        min_significance = 0.60,
        max_significance = 0.90,
        samples_significance = 0.001)
})
save(dataset.random.multiple.biopses.15.nodes, file="RData/dataset.random.multiple.biopses.15.nodes.RData")

# generate dataset for random multiple.biopses of 20 nodes
cat('dataset random multiple.biopses of 20 nodes\n')
dataset.random.multiple.biopses.20.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses(type="random",
        samples_num = sample_sizes_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate,
        nodes = 20,
        min_significance = 0.60,
        max_significance = 0.90,
        samples_significance = 0.001)
})
save(dataset.random.multiple.biopses.20.nodes, file="RData/dataset.random.multiple.biopses.20.nodes.RData")

stopCluster(cl)

print('scite multiple')
create.scite.input(dataset.random.multiple.biopses.5.nodes, 'multiple', 'random_5', scite.sd, pass.error.rates = FALSE)
create.scite.input(dataset.random.multiple.biopses.10.nodes, 'multiple', 'random_10', scite.sd, pass.error.rates = FALSE)
create.scite.input(dataset.random.multiple.biopses.15.nodes, 'multiple', 'random_15', scite.sd, pass.error.rates = FALSE)
create.scite.input(dataset.random.multiple.biopses.20.nodes, 'multiple', 'random_20', scite.sd, pass.error.rates = FALSE)
