# set the working directory
setwd("./")

# source the needed script
library("TRONCO")
source('./sample.polyclonal.trees.R')
source('./perform.polyclonal.tree.experiments.R')
source('./statistics.R')

# set the seed
set.seed(12345)

# setting of the experiments
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
sample_levels = 10
noise_levels = 8
my_algorithms = c("capri","caprese","edmonds","chowliu","prim")
my_regularizators = c("no.reg.res","loglik.res","aic.res","bic.res")

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
sampele_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
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
sampele_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)
e_pos_multiple_biopses = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
e_neg_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
low_probs_multiple_biopses = c(1.0,0.6,0.6,0.7,0.8,0.7)
medium_probs_multiple_biopses = c(1.0,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6)
high_probs_multiple_biopses = c(1.0,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6,0.8,0.7,0.6,0.7,0.7,0.8)

# run the experiments for single cells low
experiments.single.cells.low = sapply(my_experiments,function(x){run.experiments.single.cells ("low",low_true_tree,sampele_sizes_single_cells,low_probs_single_cells,e_pos_single_cells,e_neg_single_cells)})
save(experiments.single.cells.low,file="RData/experiments.single.cells.low.RData")

# statistics for the experiments for single cells low
experiments.single.cells.low.stats = get.stats(experiments.single.cells.low,my_algorithms,my_regularizators,sample_levels,noise_levels,number_experiments)
save(experiments.single.cells.low.stats,file="RData/experiments.single.cells.low.stats.RData")

# run the experiments for single cells medium
experiments.single.cells.medium = sapply(my_experiments,function(x){run.experiments.single.cells ("medium",medium_true_tree,sampele_sizes_single_cells,medium_probs_single_cells,e_pos_single_cells,e_neg_single_cells)})
save(experiments.single.cells.medium,file="RData/experiments.single.cells.medium.RData")

# statistics for the experiments for single cells medium
experiments.single.cells.medium.stats = get.stats(experiments.single.cells.medium,my_algorithms,my_regularizators,sample_levels,noise_levels,number_experiments)
save(experiments.single.cells.medium.stats,file="RData/experiments.single.cells.medium.stats.RData")

# run the experiments for single cells high
experiments.single.cells.high = sapply(my_experiments,function(x){run.experiments.single.cells ("high",high_true_tree,sampele_sizes_single_cells,high_probs_single_cells,e_pos_single_cells,e_neg_single_cells)})
save(experiments.single.cells.high,file="RData/experiments.single.cells.high.RData")

# statistics for the experiments for single cells high
experiments.single.cells.high.stats = get.stats(experiments.single.cells.high,my_algorithms,my_regularizators,sample_levels,noise_levels,number_experiments)
save(experiments.single.cells.high.stats,file="RData/experiments.single.cells.high.stats.RData")

# run the experiments for multiple biopses low
experiments.multiple.biopses.low = sapply(my_experiments,function(x){run.experiments.multiple.biopses ("low",low_true_tree,sampele_sizes_multiple_biopses,clones_num_sampling_low, low_probs_multiple_biopses,e_pos_multiple_biopses,e_neg_multiple_biopses,wild_type_rate)})
save(experiments.multiple.biopses.low,file="RData/experiments.multiple.biopses.low.RData")

# statistics for the experiments for multiple biopses low
experiments.multiple.biopses.low.stats = get.stats(experiments.multiple.biopses.low,my_algorithms,my_regularizators,sample_levels,noise_levels,number_experiments)
save(experiments.multiple.biopses.low.stats,file="RData/experiments.multiple.biopses.low.stats.RData")

# run the experiments for multiple biopses medium
experiments.multiple.biopses.medium = sapply(my_experiments,function(x){run.experiments.multiple.biopses ("medium",medium_true_tree,sampele_sizes_multiple_biopses,clones_num_sampling_medium,medium_probs_multiple_biopses,e_pos_multiple_biopses,e_neg_multiple_biopses,wild_type_rate)})
save(experiments.multiple.biopses.medium,file="RData/experiments.multiple.biopses.medium.RData")

# statistics for the experiments for multiple biopses medium
experiments.multiple.biopses.medium.stats = get.stats(experiments.multiple.biopses.medium,my_algorithms,my_regularizators,sample_levels,noise_levels,number_experiments)
save(experiments.multiple.biopses.medium.stats,file="RData/experiments.multiple.biopses.medium.stats.RData")

# run the experiments for multiple biopses high
experiments.multiple.biopses.high = sapply(my_experiments,function(x){run.experiments.multiple.biopses ("high",high_true_tree,sampele_sizes_multiple_biopses,clones_num_sampling_high, high_probs_multiple_biopses,e_pos_multiple_biopses,e_neg_multiple_biopses,wild_type_rate)})
save(experiments.multiple.biopses.high,file="RData/experiments.multiple.biopses.high.RData")

# statistics for experiments for multiple biopses high
experiments.multiple.biopses.high.stats = get.stats(experiments.multiple.biopses.high,my_algorithms,my_regularizators,sample_levels,noise_levels,number_experiments)
save(experiments.multiple.biopses.high.stats,file="RData/experiments.multiple.biopses.high.stats.RData")
