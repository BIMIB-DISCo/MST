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

if (!dir.exists('RData')) {
    dir.create('RData')
}


source('generate.dataset.R')
source('generate.scite.input.R')

# setting of the experiments
seed = 12345
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
sample_levels = 10
noise_levels = 8
cores.ratio = 1
scite.sd = 0

# set the true trees
low_true_tree = array(0,c(6,6))
low_true_tree[1,2] = 1
low_true_tree[1,3] = 1
low_true_tree[1,4] = 1
low_true_tree[2,5] = 1
low_true_tree[3,6] = 1

medium_true_tree = array(0,c(11,11))
medium_true_tree[1,2] = 1
medium_true_tree[1,3] = 1
medium_true_tree[1,4] = 1
medium_true_tree[2,5] = 1
medium_true_tree[4,6] = 1
medium_true_tree[5,7] = 1
medium_true_tree[5,8] = 1
medium_true_tree[6,9] = 1
medium_true_tree[6,10] = 1
medium_true_tree[6,11] = 1

high_true_tree = array(0,c(17,17))
high_true_tree[1,2] = 1
high_true_tree[1,3] = 1
high_true_tree[1,4] = 1
high_true_tree[2,5] = 1
high_true_tree[4,6] = 1
high_true_tree[5,7] = 1
high_true_tree[5,8] = 1
high_true_tree[6,9] = 1
high_true_tree[6,10] = 1
high_true_tree[6,11] = 1
high_true_tree[7,12] = 1
high_true_tree[9,13] = 1
high_true_tree[9,14] = 1
high_true_tree[9,15] = 1
high_true_tree[11,16] = 1
high_true_tree[11,17] = 1

# setting for single cell
sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
e_pos_single_cells = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
e_neg_single_cells = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
low_probs_single_cells = c(0.9,0.6,0.6,0.7,0.8,0.7)
medium_probs_single_cells = c(0.9,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6)
high_probs_single_cells = c(0.9,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6,0.8,0.7,0.6,0.7,0.7,0.8)

# setting for multiple biopses
clones_num_sampling_low = 3
clones_num_sampling_medium = 5
clones_num_sampling_high = 8
wild_type_rate = 0.1
sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)
e_pos_multiple_biopses = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
e_neg_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
low_probs_multiple_biopses = c(1.0,0.6,0.6,0.7,0.8,0.7)
medium_probs_multiple_biopses = c(1.0,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6)
high_probs_multiple_biopses = c(1.0,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6,0.8,0.7,0.6,0.7,0.7,0.8)


cores = as.integer(cores.ratio * (detectCores() - 1))
if (cores < 1) {
    cores = 1
}

cl = makeCluster(cores)
clusterEvalQ(cl, source('generate.dataset.R'))
clusterEvalQ(cl, source('generate.sample.polyclonal.trees.R'))
clusterExport(cl, c('sample_sizes_single_cells',
    'low_probs_single_cells',
    'medium_probs_single_cells',
    'high_probs_single_cells',
    'e_pos_single_cells',
    'e_neg_single_cells'))
clusterExport(cl, c('low_true_tree',
    'medium_true_tree',
    'high_true_tree'))
clusterExport(cl, c('clones_num_sampling_low',
    'clones_num_sampling_medium',
    'clones_num_sampling_high',
    'wild_type_rate',
    'sample_sizes_multiple_biopses',
    'e_pos_multiple_biopses',
    'e_neg_multiple_biopses',
    'low_probs_multiple_biopses',
    'medium_probs_multiple_biopses',
    'high_probs_multiple_biopses'))
clusterSetRNGStream(cl, iseed = seed)

cat('Using', cores, 'cores via "parallel" \n')

# generate dataset for single cells low
cat('dataset single cells low\n')
dataset.single.cells.low = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells("low",
        low_true_tree,
        sample_sizes_single_cells,
        low_probs_single_cells,
        e_pos_single_cells,
        e_neg_single_cells)
})
save(dataset.single.cells.low, file="RData/dataset.single.cells.low.RData")

# generate dataset for single cells medium
cat('dataset single cells medium\n')
dataset.single.cells.medium = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells("medium",
        medium_true_tree,
        sample_sizes_single_cells,
        medium_probs_single_cells,
        e_pos_single_cells,
        e_neg_single_cells)
})
save(dataset.single.cells.medium, file="RData/dataset.single.cells.medium.RData")

# generate dataset for single cells high
cat('dataset single cells high\n')
dataset.single.cells.high = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells("high",
        high_true_tree,
        sample_sizes_single_cells,
        high_probs_single_cells,
        e_pos_single_cells,
        e_neg_single_cells)
})
save(dataset.single.cells.high, file="RData/dataset.single.cells.high.RData")

# generate dataset for multiple biopses low
cat('dataset multiple biopses low\n')
dataset.multiple.biopses.low = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("low",
        low_true_tree,
        sample_sizes_multiple_biopses,
        clones_num_sampling_low,
        low_probs_multiple_biopses,
        e_pos_multiple_biopses,
        e_neg_multiple_biopses,
        wild_type_rate)
})
save(dataset.multiple.biopses.low,file="RData/dataset.multiple.biopses.low.RData")

# generate dataset for multiple biopses medium
cat('dataset multiple biopses medium\n')
dataset.multiple.biopses.medium = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("medium",
        medium_true_tree,
        sample_sizes_multiple_biopses,
        clones_num_sampling_medium,
        medium_probs_multiple_biopses,
        e_pos_multiple_biopses,
        e_neg_multiple_biopses,
        wild_type_rate)
})
save(dataset.multiple.biopses.medium,file="RData/dataset.multiple.biopses.medium.RData")

# generate dataset for multiple biopses high
cat('dataset multiple biopses high\n')
dataset.multiple.biopses.high = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("high",
        high_true_tree,
        sample_sizes_multiple_biopses,
        clones_num_sampling_high,
        high_probs_multiple_biopses,
        e_pos_multiple_biopses,
        e_neg_multiple_biopses,
        wild_type_rate)
})
save(dataset.multiple.biopses.high,file="RData/dataset.multiple.biopses.high.RData")

stopCluster(cl)

#generate scite dataset
if (dir.exists('scite_input')) {
    unlink('scite_input', recursive = TRUE)
    unlink('scite.script.*')
}

create.scite.input(dataset.single.cells.low, 'single', 'low', scite.sd)
create.scite.input(dataset.single.cells.medium, 'single', 'medium', scite.sd)
create.scite.input(dataset.single.cells.high, 'single', 'high', scite.sd)
create.scite.input(dataset.multiple.biopses.low, 'multiple', 'low', scite.sd)
create.scite.input(dataset.multiple.biopses.medium, 'multiple', 'medium', scite.sd)
create.scite.input(dataset.multiple.biopses.high, 'multiple', 'high', scite.sd)
