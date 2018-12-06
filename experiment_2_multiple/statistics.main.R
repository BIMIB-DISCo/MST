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

load('RData/experiments.multiple.biopses.low.scite.RData')
load('RData/experiments.multiple.biopses.medium.scite.RData')
load('RData/experiments.multiple.biopses.high.scite.RData')

experiments.multiple.biopses.low.scite.stats =
    get.stats(experiments.multiple.biopses.low.scite)
save(experiments.multiple.biopses.low.scite.stats,
     file = "RData/experiments.multiple.biopses.low.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.low.scite.stats, 'multiple', 'low')

experiments.multiple.biopses.medium.scite.stats =
    get.stats(experiments.multiple.biopses.medium.scite)
save(experiments.multiple.biopses.medium.scite.stats,
     file = "RData/experiments.multiple.biopses.medium.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.medium.scite.stats, 'multiple', 'medium')

experiments.multiple.biopses.high.scite.stats =
    get.stats(experiments.multiple.biopses.high.scite)
save(experiments.multiple.biopses.high.scite.stats,
     file = "RData/experiments.multiple.biopses.high.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.high.scite.stats, 'multiple', 'high')

### end of file -- statistics.main.R
