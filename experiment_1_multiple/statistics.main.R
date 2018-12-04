##############################################################################
###
### MST
###
### Statistics Main
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

load('RData/experiments.multiple.biopses.clean.scite.RData')
load('RData/experiments.multiple.biopses.convergent.scite.RData')
load('RData/experiments.multiple.biopses.random.columns.scite.RData')
load('RData/experiments.multiple.biopses.random.forest.scite.RData')

experiments.multiple.biopses.clean.scite.stats =
    get.stats(experiments.multiple.biopses.clean.scite)
save(experiments.multiple.biopses.clean.scite.stats,
     file = "RData/experiments.multiple.biopses.clean.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.clean.scite.stats,
            'multiple',
            'clean',
            c(20))

experiments.multiple.biopses.convergent.scite.stats =
    get.stats(experiments.multiple.biopses.convergent.scite)
save(experiments.multiple.biopses.convergent.scite.stats,
     file = "RData/experiments.multiple.biopses.convergent.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.convergent.scite.stats,
            'multiple',
            'convergent',
            c(20))

experiments.multiple.biopses.random.columns.scite.stats =
    get.stats(experiments.multiple.biopses.random.columns.scite)
save(experiments.multiple.biopses.random.columns.scite.stats,
     file = "RData/experiments.multiple.biopses.random.columns.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.random.columns.scite.stats,
            'multiple',
            'random_columns',
            c(20))

experiments.multiple.biopses.random.forest.scite.stats =
    get.stats(experiments.multiple.biopses.random.forest.scite)
save(experiments.multiple.biopses.random.forest.scite.stats,
     file = "RData/experiments.multiple.biopses.random.forest.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.random.forest.scite.stats,
            'multiple',
            'random_forest',
            c(20))

### end of file -- statistics.main.R
