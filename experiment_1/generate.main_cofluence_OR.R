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
source('../generate.sample.polyclonal.trees.R')

# setting of the experiments
seed = 12345
number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)

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


# setting for single cell
sample_sizes_single_cells = c(75)
e_pos_single_cells = c(0.005)
e_neg_single_cells = c(0.050)

probs_single_cells = c(0.9,0.6,0.6,0.7,0.8,0.7,0.9,0.8,0.7,0.7,0.6)

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
clusterExport(cl, c('sample_sizes_single_cells',
    'probs_single_cells',
    'e_pos_single_cells',
    'e_neg_single_cells'))
clusterExport(cl, c('convergent_true_tree'))
clusterSetRNGStream(cl, iseed = seed)


cat('Using', cores, 'cores via "parallel" \n')

# generate dataset for single cells convergent
cat('dataset single cells convergent\n')
dataset.single.cells.convergent = parSapply(cl, my_experiments, function(x){
    generate.dataset.single.cells.convergent("medium",
        convergent_true_tree,
        sample_sizes_single_cells,
        probs_single_cells,
        e_pos_single_cells,
        e_neg_single_cells)
})
dataset.single.cells.convergent = t(as.matrix(dataset.single.cells.convergent, byrow = TRUE))
colnames(dataset.single.cells.convergent) = names(my_experiments)
rownames(dataset.single.cells.convergent) = 75
save(dataset.single.cells.convergent, file="RData/dataset.single.cells.convergent_OR.RData")

stopCluster(cl)
