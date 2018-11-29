##############################################################################
###
### MST
###
### Reconstruction (slm)
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

### Source the needed script

library(TRONCO)

source('reconstruct.run.slmh.R')


### Setting of the experiments

seed = 12345
cores = 3

### Set the true trees

low_true_tree = array(0,c(6, 6))
low_true_tree[1, 2] = 1
low_true_tree[1, 3] = 1
low_true_tree[1, 4] = 1
low_true_tree[2, 5] = 1
low_true_tree[3, 6] = 1

medium_true_tree = array(0,c(11, 11))
medium_true_tree[1, 2] = 1
medium_true_tree[1, 3] = 1
medium_true_tree[1, 4] = 1
medium_true_tree[2, 5] = 1
medium_true_tree[4, 6] = 1
medium_true_tree[5, 7] = 1
medium_true_tree[5, 8] = 1
medium_true_tree[6, 9] = 1
medium_true_tree[6, 10] = 1
medium_true_tree[6, 11] = 1


### Setting for single cell

sample_sizes_single_cells = c(10, 50, 100, 250)
epos = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
eneg = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
low_probs_single_cells = c(0.9, 0.6, 0.6, 0.7, 0.8, 0.7)
medium_probs_single_cells = c(0.9, 0.6, 0.6, 0.7, 0.8, 0.7, 0.9, 0.8, 0.7, 0.7, 0.6)

load('RData/dataset.single.cells.low.final.RData')
load('RData/dataset.single.cells.medium.final.RData')


### Generate dataset for single cells low

cat('result single cells low\n')
result.single.cells.low = expand.input(dataset.single.cells.low.final, low_true_tree, seed, cores, epos, eneg)
save(result.single.cells.low, file="RData/result.single.cells.low.RData")


### Generate dataset for single cells medium

cat('result single cells medium\n')
result.single.cells.medium = expand.input(dataset.single.cells.medium.final, medium_true_tree, seed, cores, epos, eneg)
save(result.single.cells.medium, file="RData/result.single.cells.medium.RData")


### end of file -- reconstruct.main.slm.R
