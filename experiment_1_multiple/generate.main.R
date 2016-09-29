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

source('../generate.dataset.R')
source('../generate.scite.input.R')
source('../generate.sample.polyclonal.trees.R')

add.columns <-function(x, ncols) {
    dataset = x$dataset
    true_tree = x$true_tree

    nrow = nrow(dataset)
    
    sums = colSums(dataset)
    marginal = sums/nrow(dataset)
    marginal.mean = mean(marginal)
    inf.limit = max(0.1, marginal.mean - 0.15)
    sup.limit = min(0.9, marginal.mean + 0.15)

    for (i in 1:ncols) {
        # compute the marginal rate
        alpha = runif(1, inf.limit, sup.limit)
        col = rep(0, nrow)
        # add if x < alpha flip 0 -> 1
        col[which(runif(nrow) < alpha)] = 1
        dataset = cbind(dataset, col)
        colnames(dataset) = NULL
        # append an empty row and an empty columt to true_tree
        dim.true_tree = nrow(true_tree)
        true_tree = cbind(true_tree, rep(0, dim.true_tree))
        true_tree = rbind(true_tree, rep(0, dim.true_tree + 1))
    }
    x$dataset = dataset
    x$true_tree = true_tree
    return(x)
}

add.random.columns <- function(datasets, ncols) {
    for(i in 1:nrow(datasets)) {
        for (j in 1:ncol(datasets)) {
            cat((((i - 1) * ncol(datasets)) + j) , '/', nrow(datasets) * ncol(datasets), '\n')
            single.experiment = datasets[[i,j]] 
            new.dataset = lapply(single.experiment, function(x){
                add.columns(x, ncols)
            })
            for (k in 1:length(new.dataset)) {
                datasets[[i,j]][[k]] = new.dataset[[k]]
            }
        }
    }
    return(datasets)

}

# setting of the experiments
seed = 12345
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)

clean_true_tree = array(0,c(11,11))
clean_true_tree[1,2] = 1
clean_true_tree[1,3] = 1
clean_true_tree[1,4] = 1
clean_true_tree[2,5] = 1
clean_true_tree[4,6] = 1
clean_true_tree[5,7] = 1
clean_true_tree[5,8] = 1
clean_true_tree[6,9] = 1
clean_true_tree[6,10] = 1
clean_true_tree[6,11] = 1

convergent_true_tree = array(0,c(8,8))
convergent_true_tree[1,2] = 1
convergent_true_tree[1,3] = 1
convergent_true_tree[1,4] = 1
convergent_true_tree[2,5] = 1
convergent_true_tree[4,6] = 1
convergent_true_tree[3,6] = 1
convergent_true_tree[3,5] = 1
convergent_true_tree[5,7] = 1
convergent_true_tree[6,8] = 1

# setting for multiple biopses
clones_num_sampling = 5
wild_type_rate = 0.1
#sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)
sample_sizes_multiple_biopses = c(20)
#e_pos_multiple_biopses = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
#e_neg_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
e_pos_multiple_biopses = c(0.050)
e_neg_multiple_biopses = c(0.050)

probs_multiple_biopses = c(1.0,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6)

available.cores = detectCores()

if (available.cores > 8) {
    cores = 8
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}

# create the cluster

cl = makeCluster(cores)
clusterEvalQ(cl, source('../generate.dataset.R'))
clusterEvalQ(cl, source('../generate.sample.polyclonal.trees.R'))
clusterExport(cl, c('sample_sizes_multiple_biopses',
    'probs_multiple_biopses',
    'e_pos_multiple_biopses',
    'e_neg_multiple_biopses'))
clusterExport(cl, c('clones_num_sampling',
    'clean_true_tree',
    'convergent_true_tree',
    'wild_type_rate'))
clusterSetRNGStream(cl, iseed = seed)

cat('Using', cores, 'cores via "parallel" \n')

# generate dataset for multiple biopses medium
cat('dataset multiple biopses low\n')
dataset.multiple.biopses.clean = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("medium",
        clean_true_tree,
        sample_sizes_multiple_biopses,
        clones_num_sampling,
        probs_multiple_biopses,
        e_pos_multiple_biopses,
        e_neg_multiple_biopses,
        wild_type_rate)
})
dataset.multiple.biopses.clean = t(as.matrix(dataset.multiple.biopses.clean, byrow = TRUE))
colnames(dataset.multiple.biopses.clean) = names(my_experiments)
rownames(dataset.multiple.biopses.clean) = 20
save(dataset.multiple.biopses.clean, file="RData/dataset.multiple.biopses.clean.RData")

# generate dataset for single cells convergent
cat('dataset multiple biopses convergent\n')
dataset.multiple.biopses.convergent = parSapply(cl, my_experiments, function(x){
    generate.dataset.multiple.biopses("medium",
        convergent_true_tree,
        sample_sizes_multiple_biopses,
        clones_num_sampling,
        probs_multiple_biopses,
        e_pos_multiple_biopses,
        e_neg_multiple_biopses,
        wild_type_rate)
})
dataset.multiple.biopses.convergent = t(as.matrix(dataset.multiple.biopses.convergent, byrow = TRUE))
colnames(dataset.multiple.biopses.convergent) = names(my_experiments)
rownames(dataset.multiple.biopses.convergent) = 20
save(dataset.multiple.biopses.convergent, file="RData/dataset.multiple.biopses.convergent.RData")

# generate dataset for multiple biopses random column
cat('dataset multiple biopses random columns\n')
dataset.multiple.biopses.random.columns = add.random.columns(dataset.multiple.biopses.clean, 2)
save(dataset.multiple.biopses.random.columns, file = 'RData/dataset.multiple.biopses.random.columns.RData')

# generate dataset for single cells convergent
cat('dataset multiple biopses forest\n')
dataset.multiple.biopses.random.forest = sapply(my_experiments, function(x){
    generate.dataset.multiple.biopses("random_forest_fixed",
        samples_num = sample_sizes_multiple_biopses,
        e_pos = e_pos_multiple_biopses,
        e_neg = e_neg_multiple_biopses,
        wild_type = wild_type_rate)
})
dataset.multiple.biopses.random.forest = t(as.matrix(dataset.multiple.biopses.random.forest, byrow = TRUE))
colnames(dataset.multiple.biopses.random.forest) = names(my_experiments)
rownames(dataset.multiple.biopses.random.forest) = 20
save(dataset.multiple.biopses.random.forest, file="RData/dataset.multiple.biopses.random.forest.RData")

stopCluster(cl)

#generate scite dataset
if (dir.exists('scite_input')) {
    unlink('scite_input', recursive = TRUE)
    unlink('scite.script.*')
}

cat('scite clean\n')
create.scite.input(dataset.multiple.biopses.clean, 'multiple', 'clean', scite.sd, pass.error.rates = FALSE)
cat('scite convergent\n')
create.scite.input(dataset.multiple.biopses.convergent, 'multiple', 'convergent', scite.sd)
cat('scite random_columns\n')
create.scite.input(dataset.multiple.biopses.random.columns, 'multiple', 'random_columns', scite.sd)
cat('scite random_columns\n')
create.scite.input(dataset.multiple.biopses.random.forest, 'multiple', 'random_forest_fixed', scite.sd)

