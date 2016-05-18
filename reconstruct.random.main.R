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
my_algorithms = c("capri","caprese","edmonds","mle","chowliu","prim", "scite")
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
save(result.random.single.cells.5.nodes, file="RData/result.random.single.cells.5.nodes.RData")

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

#stop the cluster
stopCluster(cl)


#### please, run SCITE!!!
if (!dir.exists('scite_output')) {
    stop('run SCITE first!')
}




library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiments.random.single.cells.5.nodes.scite = import.scite.output(result.random.single.cells.5.nodes, 'single', 'random_5')
save(experiments.random.single.cells.5.nodes.scite, file = 'RData/experiments.random.single.cells.5.nodes.scite.RData')

experiments.random.single.cells.10.nodes.scite = import.scite.output(result.random.single.cells.10.nodes, 'single', 'random_10')
save(experiments.random.single.cells.10.nodes.scite, file = 'RData/experiments.random.single.cells.10.nodes.scite.RData')

experiments.random.single.cells.15.nodes.scite = import.scite.output(result.random.single.cells.15.nodes, 'single', 'random_15')
save(experiments.random.single.cells.15.nodes.scite, file = 'RData/experiments.random.single.cells.15.nodes.scite.RData')

experiments.random.single.cells.20.nodes.scite = import.scite.output(result.random.single.cells.20.nodes, 'single', 'random_20')
save(experiments.random.single.cells.20.nodes.scite, file = 'RData/experiments.random.single.cells.20.nodes.scite.RData')

experiments.random.multiple.biopses.5.nodes.scite = import.scite.output(result.random.multiple.biopses.5.nodes, 'multiple', 'random_5')
save(experiments.random.multiple.biopses.5.nodes.scite, file = 'RData/experiments.random.multiple.biopses.5.nodes.scite.RData')

experiments.random.multiple.biopses.10.nodes.scite = import.scite.output(result.random.multiple.biopses.10.nodes, 'multiple', 'random_10')
save(experiments.random.multiple.biopses.10.nodes.scite, file = 'RData/experiments.random.multiple.biopses.10.nodes.scite.RData')

experiments.random.multiple.biopses.15.nodes.scite = import.scite.output(result.random.multiple.biopses.15.nodes, 'multiple', 'random_15')
save(experiments.random.multiple.biopses.15.nodes.scite, file = 'RData/experiments.random.multiple.biopses.15.nodes.scite.RData')

experiments.random.multiple.biopses.20.nodes.scite = import.scite.output(result.random.multiple.biopses.20.nodes, 'multiple', 'random_20')
save(experiments.random.multiple.biopses.20.nodes.scite, file = 'RData/experiments.random.multiple.biopses.20.nodes.scite.RData')

