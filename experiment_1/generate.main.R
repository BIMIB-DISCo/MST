##############################################################################
###
### MST
###
### Generate Main
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

## source the needed script

library(TRONCO)
library(parallel)

if (!dir.exists('RData')) {
    dir.create('RData')
}

source('../generate.dataset.R')
source('../generate.scite.input.R')


add.columns <- function(x, ncols) {
    dataset = x$dataset
    true_tree = x$true_tree

    nrow = nrow(dataset)
    
    sums = colSums(dataset)
    marginal = sums/nrow(dataset)
    marginal.mean = mean(marginal)
    inf.limit = max(0.1, marginal.mean - 0.15)
    sup.limit = min(0.9, marginal.mean + 0.15)

    for (i in 1:ncols) {
        ## compute the marginal rate
        alpha = runif(1, inf.limit, sup.limit)
        col = rep(0, nrow)
        ## add if x < alpha flip 0 -> 1
        col[which(runif(nrow) < alpha)] = 1
        dataset = cbind(dataset, col)
        colnames(dataset) = NULL
        ## append an empty row and an empty columt to true_tree
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
            cat((((i - 1) * ncol(datasets)) + j) , '/',
                nrow(datasets) * ncol(datasets), '\n')
            single.experiment = datasets[[i, j]] 
            new.dataset = lapply(single.experiment,
                                 function(x) {
                                     add.columns(x, ncols)
                                 })
            for (k in 1:length(new.dataset)) {
                datasets[[i, j]][[k]] = new.dataset[[k]]
            }
        }
    }
    return(datasets)
}


## setting for the experiments

seed = 12345
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment", my_experiments)

clean_true_tree = array(0, c(11, 11))
clean_true_tree[1, 2] = 1
clean_true_tree[1, 3] = 1
clean_true_tree[1, 4] = 1
clean_true_tree[2, 5] = 1
clean_true_tree[4, 6] = 1
clean_true_tree[5, 7] = 1
clean_true_tree[5, 8] = 1
clean_true_tree[6, 9] = 1
clean_true_tree[6, 10] = 1
clean_true_tree[6, 11] = 1


convergent_true_tree = array(0, c(8, 8))
convergent_true_tree[1, 2] = 1
convergent_true_tree[1, 3] = 1
convergent_true_tree[1, 4] = 1
convergent_true_tree[2, 5] = 1
convergent_true_tree[4, 6] = 1
convergent_true_tree[3, 6] = 1
convergent_true_tree[3, 5] = 1
convergent_true_tree[5, 7] = 1
convergent_true_tree[6, 8] = 1


## Setting for single cell

sample_sizes_single_cells = c(75)
##e_pos_single_cells = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
##e_neg_single_cells = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
e_pos_single_cells = c(0.005)
e_neg_single_cells = c(0.050)



probs_single_cells = c(0.9, 0.6, 0.6, 0.7, 0.8, 0.7, 0.9, 0.8, 0.7, 0.7, 0.6)

available.cores = detectCores()

if (available.cores > 8) {
    cores = 8
} else if (available.cores > 1) {
    cores = available.cores - 1
} else {
    cores = 1
}

## create the cluster

cl = makeCluster(cores)
clusterEvalQ(cl, source('../generate.dataset.R'))
clusterEvalQ(cl, source('../generate.sample.polyclonal.trees.R'))
clusterExport(cl, c('sample_sizes_single_cells', 
                    'probs_single_cells', 
                    'e_pos_single_cells', 
                    'e_neg_single_cells'))
clusterExport(cl, c('clean_true_tree', 
    'convergent_true_tree'))
clusterSetRNGStream(cl, iseed = seed)


cat('Using', cores, 'cores via "parallel" \n')


## Generate dataset for single cells clean

cat('dataset single cells clean\n')
dataset.single.cells.clean =
    parSapply(cl,
              my_experiments,
              function(x) {
                  generate.dataset.single.cells("medium", 
                                                clean_true_tree, 
                                                sample_sizes_single_cells, 
                                                probs_single_cells, 
                                                e_pos_single_cells, 
                                                e_neg_single_cells)
              })
dataset.single.cells.clean = t(as.matrix(dataset.single.cells.clean, byrow = TRUE))
colnames(dataset.single.cells.clean) = names(my_experiments)
rownames(dataset.single.cells.clean) = 75
save(dataset.single.cells.clean, file = "RData/dataset.single.cells.clean.RData")


## Generate dataset for single cells convergent

cat('dataset single cells convergent\n')
dataset.single.cells.convergent =
    parSapply(cl,
              my_experiments,
              function(x) {
                  generate.dataset.single.cells("medium", 
                                                convergent_true_tree, 
                                                sample_sizes_single_cells, 
                                                probs_single_cells, 
                                                e_pos_single_cells, 
                                                e_neg_single_cells)
              })
dataset.single.cells.convergent = t(as.matrix(dataset.single.cells.convergent, byrow = TRUE))
colnames(dataset.single.cells.convergent) = names(my_experiments)
rownames(dataset.single.cells.convergent) = 75
save(dataset.single.cells.convergent, file = "RData/dataset.single.cells.convergent.RData")


## Generate dataset for single cells random column

cat('dataset single cells random columns\n')
dataset.single.cells.random.columns = add.random.columns(dataset.single.cells.clean, 2)
save(dataset.single.cells.random.columns, file = 'RData/dataset.single.cells.random.columns.RData')


## Generate dataset for random single cells of 5 nodes

cat('dataset random single cells forest of 20 nodes\n')
dataset.single.cells.random.forest =
    sapply(my_experiments,
           function(x) {
               generate.dataset.single.cells(type = "random_forest_fixed", 
                                             samples_num = sample_sizes_single_cells, 
                                             e_pos = e_pos_single_cells, 
                                             e_neg = e_neg_single_cells, 
                                             min_significance = 0.70, 
                                             max_significance = 0.90, 
                                             samples_significance = 0.001)
           })
dataset.single.cells.random.forest = t(as.matrix(dataset.single.cells.random.forest, byrow = TRUE))
colnames(dataset.single.cells.random.forest) = names(my_experiments)
rownames(dataset.single.cells.random.forest) = 75
save(dataset.single.cells.random.forest, file = "RData/dataset.single.cells.random.forest.RData")


stopCluster(cl)


## Generate scite dataset
if (dir.exists('scite_input')) {
    unlink('scite_input', recursive = TRUE)
    unlink('scite.script.*')
}

cat('scite clean\n')
create.scite.input(dataset.single.cells.clean, 'single', 'clean', scite.sd)
cat('scite convergent\n')
create.scite.input(dataset.single.cells.convergent, 'single', 'convergent', scite.sd)
cat('scite random_columns\n')
create.scite.input(dataset.single.cells.random.columns, 'single', 'random_columns', scite.sd)
cat('scite random_forest\n')
create.scite.input(dataset.single.cells.random.forest, 'single', 'random_forest_fixed', scite.sd)


### end of file -- generate.main.R
