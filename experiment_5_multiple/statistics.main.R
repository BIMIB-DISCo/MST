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


check.fallback <- function(dataset) {
    samples = c(5, 7, 10, 20, 50)
    error = c(0.000, 0.05, 0.1, 0.15, 0.20)
    dataset = dataset$fallback_edmonds$gabow$pmi.no.reg
    for (i in 1:nrow(dataset)) {
        for(j in 1:ncol(dataset)) {
            obj = dataset[[i,j]]
            if (obj$sum > 0) {
                set = which(obj$values == 1)
                cat('sample: ', samples[[i]],
                    ' - noise: ', error[[j]],
                    ' - sum(fallback_edmonds): ', obj$sum,
                    ' - exp: ', set, '\n')
            }
        }
    }
}


check.fallback(experiments.multiple.biopses.low.scite.stats)
check.fallback(experiments.multiple.biopses.medium.scite.stats)
check.fallback(experiments.multiple.biopses.high.scite.stats)

### end of file -- statistics.main.R
