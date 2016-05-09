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


source('generate.dataset.R')
source('generate.scite.input.R')

# setting of the experiments
seed = 12345
number_experiments = 1#100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
sample_levels = 10
noise_levels = 8
cores.ratio = 1
scite.sd = 0

# setting for single cell
sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
e_pos_single_cells = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
e_neg_single_cells = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)

# setting for multiple biopses
wild_type_rate = 0.1
sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)
e_pos_multiple_biopses = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
e_neg_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)


cores = as.integer(cores.ratio * (detectCores() - 1))
if (cores < 1) {
    cores = 1
}

cl = makeCluster(cores)
clusterEvalQ(cl, source('generate.dataset.R'))
clusterEvalQ(cl, source('generate.sample.polyclonal.trees.R'))
clusterExport(cl, c('sample_sizes_single_cells',
    'e_pos_single_cells',
    'e_neg_single_cells'))
clusterExport(cl, c('wild_type_rate',
    'sample_sizes_multiple_biopses',
    'e_pos_multiple_biopses',
    'e_neg_multiple_biopses'))
clusterSetRNGStream(cl, iseed = seed)

cat('Using', cores, 'cores via "parallel" \n')

# generate dataset for random single cells of 5 nodes
cat('dataset random single cells of 5 nodes\n')
dataset.random.single.cells.5.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells(type="random",
        true_tree = NULL,
        samples_num = sample_sizes_single_cells,
        low_probs_single_cells,
        e_pos = e_pos_single_cells,
        e_neg = e_neg_single_cells,
        nodes = 5,
        significance = 0.10,
        samples_significance = 0.001)
})
save(dataset.random.single.cells.5.nodes, file="RData/dataset.random.single.cells.5.nodes.RData")

# generate dataset for random single cells of 10 nodes
cat('dataset random single cells of 10 nodes\n')
dataset.random.single.cells.10.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells(type="random",
        true_tree = NULL,
        samples_num = sample_sizes_single_cells,
        low_probs_single_cells,
        e_pos = e_pos_single_cells,
        e_neg = e_neg_single_cells,
        nodes = 10,
        significance = 0.10,
        samples_significance = 0.001)
})
save(dataset.random.single.cells.10.nodes, file="RData/dataset.random.single.cells.10.nodes.RData")

# generate dataset for random single cells of 15 nodes
cat('dataset random single cells of 15 nodes\n')
dataset.random.single.cells.15.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells(type="random",
        true_tree = NULL,
        samples_num = sample_sizes_single_cells,
        low_probs_single_cells,
        e_pos = e_pos_single_cells,
        e_neg = e_neg_single_cells,
        nodes = 15,
        significance = 0.10,
        samples_significance = 0.001)
})
save(dataset.random.single.cells.15.nodes, file="RData/dataset.random.single.cells.15.nodes.RData")

# generate dataset for multiple biopses of 5 nodes
cat('dataset multiple biopses of 5 nodes\n')
dataset.multiple.biopses.5.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("random",
        true_tree = NULL,
        samples_num = sample_sizes_multiple_biopses,
        clones_num_sampling_low,
        low_probs_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate,
        nodes = 5,
        significance = 0.10,
        samples_significance = 0.001)
})
save(dataset.multiple.biopses.5.nodes,file="RData/dataset.multiple.biopses.5.nodes.RData")

# generate dataset for multiple biopses of 10 nodes
cat('dataset multiple biopses of 10 nodes\n')
dataset.multiple.biopses.10.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("random",
        true_tree = NULL,
        samples_num = sample_sizes_multiple_biopses,
        clones_num_sampling_low,
        low_probs_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate,
        nodes = 10,
        significance = 0.10,
        samples_significance = 0.001)
})
save(dataset.multiple.biopses.10.nodes,file="RData/dataset.multiple.biopses.10.nodes.RData")

# generate dataset for multiple biopses of 15 nodes
cat('dataset multiple biopses of 15 nodes\n')
dataset.multiple.biopses.15.nodes = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("random",
        true_tree = NULL,
        samples_num = sample_sizes_multiple_biopses,
        clones_num_sampling_low,
        low_probs_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate,
        nodes = 15,
        significance = 0.10,
        samples_significance = 0.001)
})
save(dataset.multiple.biopses.15.nodes,file="RData/dataset.multiple.biopses.15.nodes.RData")


stopCluster(cl)

#generate scite dataset
if (dir.exists('scite_input')) {
    unlink('scite_input', recursive = TRUE)
    unlink('scite.script.*')
}

create.scite.input(dataset.random.single.cells.5.nodes, 'single', 'random_5', scite.sd)
create.scite.input(dataset.random.single.cells.10.nodes, 'single', 'random_10', scite.sd)
create.scite.input(dataset.random.single.cells.15.nodes, 'single', 'random_15', scite.sd)
create.scite.input(dataset.random.multiple.biopses.5.nodes, 'multiple', 'random_5', scite.sd)
create.scite.input(dataset.random.multiple.biopses.10.nodes, 'multiple', 'random_10', scite.sd)
create.scite.input(dataset.random.multiple.biopses.15.nodes, 'multiple', 'random_15', scite.sd)
