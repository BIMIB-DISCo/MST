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

source('../statistics.plot.R')
source('../statistics.compute.R')
source('../giulio.plot.R')

load('RData/experiments.random.multiple.biopses.5.nodes.scite.RData')
load('RData/experiments.random.multiple.biopses.10.nodes.scite.RData')
load('RData/experiments.random.multiple.biopses.15.nodes.scite.RData')
load('RData/experiments.random.multiple.biopses.20.nodes.scite.RData')

experiments.random.multiple.biopses.5.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.5.nodes.scite)
save(experiments.random.multiple.biopses.5.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.5.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.5.nodes.scite.stats, 'multiple', 'random_5')

experiments.random.multiple.biopses.10.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.10.nodes.scite)
save(experiments.random.multiple.biopses.10.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.10.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.10.nodes.scite.stats, 'multiple', 'random_10')

experiments.random.multiple.biopses.15.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.15.nodes.scite)
save(experiments.random.multiple.biopses.15.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.15.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.15.nodes.scite.stats, 'multiple', 'random_15')

experiments.random.multiple.biopses.20.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.20.nodes.scite)
save(experiments.random.multiple.biopses.20.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.20.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.20.nodes.scite.stats, 'multiple', 'random_20')


check.fallback = function(dataset) {
	samples = c(5, 7, 10, 20, 50)
	error = c(0.000, 0.05, 0.1, 0.15, 0.20)
	dataset = dataset$fallback_edmonds$gabow$pmi.no.reg
	for (i in 1:nrow(dataset)) {
		for(j in 1:ncol(dataset)) {
			obj = dataset[[i,j]]
			if (obj$sum > 0) {
				set = which(obj$values == 1)
				cat('sample: ', samples[[i]], ' - noise: ', error[[j]], ' - sum(fallback_edmonds): ', obj$sum, ' - exp: ', set, '\n')
			}
		}
	}
}
