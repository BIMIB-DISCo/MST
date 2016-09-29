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

source('../generate.scite.input.R')
scite.sd = 0

load('RData/dataset.multiple.biopses.low.RData')
load('RData/dataset.multiple.biopses.medium.RData')
load('RData/dataset.multiple.biopses.high.RData')

ncol.low = 1
ncol.medium = 2
ncol.high = 4

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

cat('low\n')
dataset.multiple.biopses.random.columns.low = add.random.columns(dataset.multiple.biopses.low, ncol.low)
save(dataset.multiple.biopses.random.columns.low, file = 'RData/dataset.multiple.biopses.random.columns.low.RData')
cat('medium\n')
dataset.multiple.biopses.random.columns.medium = add.random.columns(dataset.multiple.biopses.medium, ncol.medium)
save(dataset.multiple.biopses.random.columns.medium, file = 'RData/dataset.multiple.biopses.random.columns.medium.RData')
cat('high\n')
dataset.multiple.biopses.random.columns.high = add.random.columns(dataset.multiple.biopses.high, ncol.high)
save(dataset.multiple.biopses.random.columns.high, file = 'RData/dataset.multiple.biopses.random.columns.high.RData')

cat('scite low\n')
create.scite.input(dataset.multiple.biopses.random.columns.low, 'multiple', 'low', scite.sd, pass.error.rates = FALSE)
cat('scite medium\n')
create.scite.input(dataset.multiple.biopses.random.columns.medium, 'multiple', 'medium', scite.sd, pass.error.rates = FALSE)
cat('scite high\n')
create.scite.input(dataset.multiple.biopses.random.columns.high, 'multiple', 'high', scite.sd, pass.error.rates = FALSE)
