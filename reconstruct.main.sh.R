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

source('reconstruct.run.slmh.R')

# setting of the experiments
seed = 12345
cores = 3

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
sample_sizes_single_cells = c(10, 50, 100, 250)
epos = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
eneg = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
high_probs_single_cells = c(0.9,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6,0.8,0.7,0.6,0.7,0.7,0.8)


load('RData/dataset.single.cells.high.final.RData')

# generate dataset for single cells high
cat('result single cells high\n')
result.single.cells.high = expand.input(dataset.single.cells.high.final, high_true_tree, seed, cores, epos, eneg)
save(result.single.cells.high, file="RData/result.single.cells.high.RData")
