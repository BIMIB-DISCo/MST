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

source('reconstruct.run.R')
source('reconstruct.scite.import.R')

# setting of the experiments
seed = 12345
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
my_algorithms = c("capri","caprese","edmonds","chowliu","prim", "scite")
my_regularizators = c("no.reg.res","loglik.res","aic.res","bic.res")
cores.ratio = 1

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


load('RData/dataset.single.cells.low.RData')
load('RData/dataset.single.cells.medium.RData')
load('RData/dataset.single.cells.high.RData')
load('RData/dataset.multiple.biopses.low.RData')
load('RData/dataset.multiple.biopses.medium.RData')
load('RData/dataset.multiple.biopses.high.RData')

# generate dataset for single cells low
cat('result single cells low\n')
result.single.cells.low = expand.input(dataset.single.cells.low, low_true_tree, seed, cores)
save(result.single.cells.low, file="RData/result.single.cells.low.RData")

# generate dataset for single cells medium
cat('result single cells medium\n')
result.single.cells.medium = expand.input(dataset.single.cells.medium, medium_true_tree, seed, cores)
save(result.single.cells.medium, file="RData/result.single.cells.medium.RData")

# generate dataset for single cells high
cat('result single cells high\n')
result.single.cells.high = expand.input(dataset.single.cells.high, high_true_tree, seed, cores)
save(result.single.cells.high, file="RData/result.single.cells.high.RData")

# generate dataset for multiple biopses low
cat('result multiple biopses low\n')
result.multiple.biopses.low = expand.input(dataset.multiple.biopses.low, low_true_tree, seed, cores)
save(result.multiple.biopses.low,file="RData/result.multiple.biopses.low.RData")

# generate dataset for multiple biopses medium
cat('result multiple biopses medium\n')
result.multiple.biopses.medium = expand.input(dataset.multiple.biopses.medium, medium_true_tree, seed, cores)
save(result.multiple.biopses.medium,file="RData/result.multiple.biopses.medium.RData")

# generate dataset for multiple biopses high
cat('result multiple biopses high\n')
result.multiple.biopses.high = expand.input(dataset.multiple.biopses.high, high_true_tree, seed, cores)
save(result.multiple.biopses.high,file="RData/result.multiple.biopses.high.RData")



#### please, run SCITE!!!
if (!dir.exists('scite_output')) {
    stop('run SCITE first!')
}

library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiments.single.cells.low.scite = import.scite.output(result.single.cells.low, 'single', 'low', low_true_tree)
save(experiments.single.cells.low.scite, file = 'RData/experiments.single.cells.low.scite.RData')

experiments.single.cells.medium.scite = import.scite.output(result.single.cells.medium, 'single', 'medium', medium_true_tree)
save(experiments.single.cells.medium.scite, file = 'RData/experiments.single.cells.medium.scite.RData')

experiments.single.cells.high.scite = import.scite.output(result.single.cells.high, 'single', 'high', high_true_tree)
save(experiments.single.cells.high.scite, file = 'RData/experiments.single.cells.high.scite.RData')

experiments.multiple.biopses.low.scite = import.scite.output(result.multiple.biopses.low, 'multiple', 'low', low_true_tree)
save(experiments.multiple.biopses.low.scite, file = 'RData/experiments.multiple.biopses.low.scite.RData')

experiments.multiple.biopses.medium.scite = import.scite.output(result.multiple.biopses.medium, 'multiple', 'medium', medium_true_tree)
save(experiments.multiple.biopses.medium.scite, file = 'RData/experiments.multiple.biopses.medium.scite.RData')

experiments.multiple.biopses.high.scite = import.scite.output(result.multiple.biopses.high, 'multiple', 'high', high_true_tree)
save(experiments.multiple.biopses.high.scite, file = 'RData/experiments.multiple.biopses.high.scite.RData')










# statistics for the experiments for single cells low
experiments.single.cells.low.stats = get.stats(experiments.single.cells.low,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.single.cells.low.stats, file="RData/experiments.single.cells.low.stats.RData")

# statistics for the experiments for single cells medium
experiments.single.cells.medium.stats = get.stats(experiments.single.cells.medium,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.single.cells.medium.stats, file="RData/experiments.single.cells.medium.stats.RData")

# statistics for the experiments for single cells high
experiments.single.cells.high.stats = get.stats(experiments.single.cells.high,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.single.cells.high.stats,file="RData/experiments.single.cells.high.stats.RData")

# statistics for the experiments for multiple biopses low
experiments.multiple.biopses.low.stats = get.stats(experiments.multiple.biopses.low,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.multiple.biopses.low.stats, file="RData/experiments.multiple.biopses.low.stats.RData")

# statistics for the experiments for multiple biopses medium
experiments.multiple.biopses.medium.stats = get.stats(experiments.multiple.biopses.medium,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.multiple.biopses.medium.stats,file="RData/experiments.multiple.biopses.medium.stats.RData")

# statistics for experiments for multiple biopses high
experiments.multiple.biopses.high.stats = get.stats(experiments.multiple.biopses.high,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.multiple.biopses.high.stats,file="RData/experiments.multiple.biopses.high.stats.RData")

#stop the cluster
stopCluster(cl)


# performance plot

performance.plot(experiments.single.cells.low.stats, 'single', 'low')
performance.plot(experiments.single.cells.medium.stats, 'single', 'medium')
performance.plot(experiments.single.cells.high.stats, 'single', 'high')
performance.plot(experiments.multiple.biopses.low.stats, 'multiple', 'low')
performance.plot(experiments.multiple.biopses.medium.stats, 'multiple', 'medium')
performance.plot(experiments.multiple.biopses.high.stats, 'multiple', 'high')

compare.performance.plot(experiments.single.cells.low.stats, 'single', 'low')
compare.performance.plot(experiments.single.cells.medium.stats, 'single', 'medium')
compare.performance.plot(experiments.single.cells.high.stats, 'single', 'high')
compare.performance.plot(experiments.multiple.biopses.low.stats, 'multiple', 'low')
compare.performance.plot(experiments.multiple.biopses.medium.stats, 'multiple', 'medium')
compare.performance.plot(experiments.multiple.biopses.high.stats, 'multiple', 'high')

compare.performance.plot.2d(experiments.single.cells.low.stats, 'single', 'low')
compare.performance.plot.2d(experiments.single.cells.medium.stats, 'single', 'medium')
compare.performance.plot.2d(experiments.single.cells.high.stats, 'single', 'high')
compare.performance.plot.2d(experiments.multiple.biopses.low.stats, 'multiple', 'low')
compare.performance.plot.2d(experiments.multiple.biopses.medium.stats, 'multiple', 'medium')
compare.performance.plot.2d(experiments.multiple.biopses.high.stats, 'multiple', 'high')


