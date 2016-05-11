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


load('RData/dataset.random.single.cells.5.nodes.RData')
load('RData/dataset.random.single.cells.10.nodes.RData')
load('RData/dataset.random.single.cells.15.nodes.RData')
load('RData/dataset.random.single.cells.20.nodes.RData')
load('RData/dataset.random.multiple.biopses.5.nodes.RData')
load('RData/dataset.random.multiple.biopses.10.nodes.RData')
load('RData/dataset.random.multiple.biopses.15.nodes.RData')
load('RData/dataset.random.multiple.biopses.20.nodes.RData')

# generate dataset for single cells 5
cat('result single cells 5\n')
result.random.single.cells.5.nodes = expand.random.input(dataset.random.single.cells.5.nodes, seed, cores)
save(result.single.cells.low, file="RData/result.random.single.cells.5.nodes.RData")

# generate dataset for single cells low
cat('result single cells 10\n')
result.random.single.cells.10.nodes = expand.random.input(dataset.random.single.cells.10.nodes, seed, cores)
save(result.random.single.cells.10.nodes, file="RData/result.random.single.cells.10.nodes.RData")

# generate dataset for single cells low
cat('result single cells 15\n')
result.random.single.cells.15.nodes = expand.random.input(dataset.random.single.cells.15.nodes, seed, cores)
save(result.random.single.cells.15.nodes, file="RData/result.random.single.cells.15.nodes.RData")

# generate dataset for single cells low
cat('result single cells 20\n')
result.random.single.cells.20.nodes = expand.random.input(dataset.random.single.cells.20.nodes, seed, cores)
save(result.random.single.cells.20.nodes, file="RData/result.random.single.cells.20.nodes.RData")

# generate dataset for single cells low
cat('result multiple biopses 5\n')
result.random.multiple.biopses.5.nodes = expand.random.input(dataset.random.multiple.biopses.5.nodes , seed, cores)
save(result.random.multiple.biopses.5.nodes, file="RData/result.random.multiple.biopses.5.nodes.RData")

# generate dataset for single cells low
cat('result multiple biopses 10\n')
result.random.multiple.biopses.10.nodes = expand.random.input(dataset.random.multiple.biopses.10.nodes , seed, cores)
save(result.random.multiple.biopses.10.nodes, file="RData/result.random.multiple.biopses.10.nodes.RData")

# generate dataset for single cells low
cat('result multiple biopses 15\n')
result.random.multiple.biopses.15.nodes = expand.random.input(dataset.random.multiple.biopses.15.nodes , seed, cores)
save(result.random.multiple.biopses.15.nodes, file="RData/result.random.multiple.biopses.15.nodes.RData")

# generate dataset for single cells low
cat('result multiple biopses 20\n')
result.random.multiple.biopses.20.nodes = expand.random.input(dataset.random.multiple.biopses.20.nodes , seed, cores)
save(result.random.multiple.biopses.20.nodes, file="RData/result.random.multiple.biopses.20.nodes.RData")


#### please, run SCITE!!!
if (!dir.exists('scite_output')) {
    stop('run SCITE first!')
}




library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiments.single.cells.low.scite = import.scite.output(result.single.cells.low, 'single', 'low')
save(experiments.single.cells.low.scite, file = 'RData/experiments.single.cells.low.scite.RData')

experiments.single.cells.medium.scite = import.scite.output(result.single.cells.medium, 'single', 'medium')
save(experiments.single.cells.medium.scite, file = 'RData/experiments.single.cells.medium.scite.RData')

experiments.single.cells.high.scite = import.scite.output(result.single.cells.high, 'single', 'high')
save(experiments.single.cells.high.scite, file = 'RData/experiments.single.cells.high.scite.RData')

experiments.multiple.biopses.low.scite = import.scite.output(result.multiple.biopses.low, 'multiple', 'low')
save(experiments.multiple.biopses.low.scite, file = 'RData/experiments.multiple.biopses.low.scite.RData')

experiments.multiple.biopses.medium.scite = import.scite.output(result.multiple.biopses.medium, 'multiple', 'medium')
save(experiments.multiple.biopses.medium.scite, file = 'RData/experiments.multiple.biopses.medium.scite.RData')

experiments.multiple.biopses.high.scite = import.scite.output(result.multiple.biopses.high, 'multiple', 'high')
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


