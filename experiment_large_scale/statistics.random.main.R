##############################################################################
###
### MST
###
### Statistics Random Main 20
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

source('../statistics.plot.R')
source('../statistics.compute.R')
source('../giulio.plot.R')

load('RData/result.random.single.cells.20.nodes.RData')

experiments.random.single.cells.20.nodes.stats =
    get.stats(result.random.single.cells.20.nodes)
save(experiments.random.single.cells.20.nodes.stats,
     file = "RData/experiments.random.single.cells.20.nodes.stats.RData")
giulio.plot(experiments.random.single.cells.20.nodes.stats,
            'single',
            'random_20',
            sample_size = c(100, 500, 1000, 5000, 10000),
            statistics = c('hamming_distance',
                           'accuracy',
                           'sensitivity',
                           'specificity',
                           'elasped_time'))


### end of file -- statistics.random.main.R
