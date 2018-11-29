##############################################################################
###
### MST
###
### Statistics Random Main
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


source('statistics.plot.R')
source('statistics.compute.R')

number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment", my_experiments)
## my_algorithms = c("capri", "caprese", "edmonds", "mle", "chowliu", "prim", "scite")
my_algorithms = c("edmonds", "mle", "scite", 'mltree', 'gabow', 'gabow_no_rising')
## my_regularizators = c("no.reg.res", "loglik.res", "aic.res", "bic.res")
my_regularizators = c("mi.no.reg.res",
                      "pmi.no.reg.res",
                      "cpmi.no.reg.res",
                      "entropy.no.reg.res",
                      "no.raising.no.reg.res",
                      "no.raising.mi.no.reg.res",
                      "no.raising.pmi.no.reg.res",
                      "no.raising.cpmi.no.reg.res",
                      "entropy.no.reg.res",
                      'no.reg.res')
sample_levels = 4
noise_levels = 8

load('RData/experiments.single.cells.low.scite.RData')
load('RData/experiments.single.cells.medium.scite.RData')
load('RData/experiments.single.cells.high.scite.RData')
load('RData/experiments.multiple.biopses.low.scite.RData')
load('RData/experiments.multiple.biopses.medium.scite.RData')
load('RData/experiments.multiple.biopses.high.scite.RData')


experiments.random.single.cells.5.nodes.scite.stats =
    get.stats(experiments.random.single.cells.5.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.single.cells.5.nodes.scite.stats,
     file = 'RData/experiments.random.single.cells.5.nodes.scite.stats.Rdata')

experiments.random.single.cells.10.nodes.scite.stats =
    get.stats(experiments.random.single.cells.10.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.single.cells.10.nodes.scite.stats,
     file = 'RData/experiments.random.single.cells.10.nodes.scite.stats.Rdata')


experiments.random.single.cells.15.nodes.scite.stats =
    get.stats(experiments.random.single.cells.15.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.single.cells.15.nodes.scite.stats,
     file = 'RData/experiments.random.single.cells.15.nodes.scite.stats.Rdata')


experiments.random.single.cells.20.nodes.scite.stats =
    get.stats(experiments.random.single.cells.20.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.single.cells.20.nodes.scite.stats,
     file = 'RData/experiments.random.single.cells.20.nodes.scite.stats.Rdata')


experiments.random.multiple.biopses.5.nodes.scite.stats =
    get.stats(experiments.random.multiple.biopses.5.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.multiple.biopses.5.nodes.scite.stats,
     file = 'RData/experiments.random.multiple.biopses.5.nodes.scite.stats.RData')


experiments.random.multiple.biopses.10.nodes.scite.stats =
    get.stats(experiments.random.multiple.biopses.10.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.multiple.biopses.10.nodes.scite.stats,
     file = 'RData/experiments.random.multiple.biopses.10.nodes.scite.stats.RData')


experiments.random.multiple.biopses.15.nodes.scite.stats =
    get.stats(experiments.random.multiple.biopses.15.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.multiple.biopses.15.nodes.scite.stats,
     file = 'RData/experiments.random.multiple.biopses.15.nodes.scite.stats.RData')


experiments.random.multiple.biopses.20.nodes.scite.stats =
    get.stats(experiments.random.multiple.biopses.20.nodes.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.random.multiple.biopses.20.nodes.scite.stats,
     file = 'RData/experiments.random.multiple.biopses.20.nodes.scite.stats.RData')


performance.plot(experiments.random.single.cells.5.nodes.scite.stats, 'single', 'random_5')
performance.plot(experiments.random.single.cells.10.nodes.scite.stats, 'single', 'random_10')
performance.plot(experiments.random.single.cells.15.nodes.scite.stats, 'single', 'random_15')
performance.plot(experiments.random.single.cells.20.nodes.scite.stats, 'single', 'random_20')

performance.plot(experiments.random.multiple.biopses.5.nodes.scite.stats, 'multiple', 'random_5')
performance.plot(experiments.random.multiple.biopses.10.nodes.scite.stats, 'multiple', 'random_10')
performance.plot(experiments.random.multiple.biopses.15.nodes.scite.stats, 'multiple', 'random_15')
performance.plot(experiments.random.multiple.biopses.20.nodes.scite.stats, 'multiple', 'random_20')

compare.performance.plot(experiments.random.single.cells.5.nodes.scite.stats, 'single', 'random_5')
compare.performance.plot(experiments.random.single.cells.10.nodes.scite.stats, 'single', 'random_10')
compare.performance.plot(experiments.random.single.cells.15.nodes.scite.stats, 'single', 'random_15')
compare.performance.plot(experiments.random.single.cells.20.nodes.scite.stats, 'single', 'random_20')

compare.performance.plot(experiments.random.multiple.biopses.5.nodes.scite.stats, 'multiple', 'random_5')
compare.performance.plot(experiments.random.multiple.biopses.10.nodes.scite.stats, 'multiple', 'random_10')
compare.performance.plot(experiments.random.multiple.biopses.15.nodes.scite.stats, 'multiple', 'random_15')
compare.performance.plot(experiments.random.multiple.biopses.20.nodes.scite.stats, 'multiple', 'random_20')

compare.performance.plot.2d(experiments.random.single.cells.5.nodes.scite.stats, 'single', 'random_5')
compare.performance.plot.2d(experiments.random.single.cells.10.nodes.scite.stats, 'single', 'random_10')
compare.performance.plot.2d(experiments.random.single.cells.15.nodes.scite.stats, 'single', 'random_15')
compare.performance.plot.2d(experiments.random.single.cells.20.nodes.scite.stats, 'single', 'random_20')

compare.performance.plot.2d(experiments.random.multiple.biopses.5.nodes.scite.stats, 'multiple', 'random_5')
compare.performance.plot.2d(experiments.random.multiple.biopses.10.nodes.scite.stats, 'multiple', 'random_10')
compare.performance.plot.2d(experiments.random.multiple.biopses.15.nodes.scite.stats, 'multiple', 'random_15')
compare.performance.plot.2d(experiments.random.multiple.biopses.20.nodes.scite.stats, 'multiple', 'random_20')



## polyclonal

experiments.single.cells.low.scite.stats =
    get.stats(experiments.single.cells.low.scite, 
              my_algorithms, 
              my_regularizators, 
              sample_levels, 
              noise_levels, 
              number_experiments)
save(experiments.single.cells.low.scite.stats,
     file = 'RData/experiments.single.cells.low.scite.stats.Rdata')



compare.performance.plot.2d(experiments.single.cells.low.scite.stats, 'single', 'low')

stop('... all done ...')

## These are to reconstruct TRONCO only

experiments.random.single.cells.5.nodes.stats =
    get.stats(result.random.single.cells.5.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)

experiments.random.single.cells.10.nodes.stats =
    get.stats(result.random.single.cells.10.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)

experiments.random.single.cells.15.nodes.stats =
    get.stats(result.random.single.cells.15.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)

experiments.random.single.cells.20.nodes.stats =
    get.stats(result.random.single.cells.20.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)

experiments.random.multiple.biopses.5.nodes.stats =
    get.stats(result.random.multiple.biopses.5.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)

experiments.random.multiple.biopses.10.nodes.stats =
    get.stats(result.random.multiple.biopses.10.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)

experiments.random.multiple.biopses.15.nodes.stats =
    get.stats(result.random.multiple.biopses.15.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)

experiments.random.multiple.biopses.20.nodes.stats =
    get.stats(result.random.multiple.biopses.20.nodes, 
              my_algorithms,
              my_regularizators, 
              sample_levels,
              noise_levels,
              number_experiments)


performance.plot(experiments.random.single.cells.5.nodes.stats, 'single', 'random_5')
performance.plot(experiments.random.single.cells.10.nodes.stats, 'single', 'random_10')
performance.plot(experiments.random.single.cells.15.nodes.stats, 'single', 'random_15')
performance.plot(experiments.random.single.cells.20.nodes.stats, 'single', 'random_20')

performance.plot(experiments.random.multiple.biopses.5.nodes.stats, 'multiple', 'random_5')
performance.plot(experiments.random.multiple.biopses.10.nodes.stats, 'multiple', 'random_10')
performance.plot(experiments.random.multiple.biopses.15.nodes.stats, 'multiple', 'random_15')
performance.plot(experiments.random.multiple.biopses.20.nodes.stats, 'multiple', 'random_20')

compare.performance.plot(experiments.random.single.cells.5.nodes.stats, 'single', 'random_5')
compare.performance.plot(experiments.random.single.cells.10.nodes.stats, 'single', 'random_10')
compare.performance.plot(experiments.random.single.cells.15.nodes.stats, 'single', 'random_15')
compare.performance.plot(experiments.random.single.cells.20.nodes.stats, 'single', 'random_20')

compare.performance.plot(experiments.random.multiple.biopses.5.nodes.stats, 'multiple', 'random_5')
compare.performance.plot(experiments.random.multiple.biopses.10.nodes.stats, 'multiple', 'random_10')
compare.performance.plot(experiments.random.multiple.biopses.15.nodes.stats, 'multiple', 'random_15')
compare.performance.plot(experiments.random.multiple.biopses.20.nodes.stats, 'multiple', 'random_20')

compare.performance.plot.2d(experiments.random.single.cells.5.nodes.stats, 'single', 'random_5')
compare.performance.plot.2d(experiments.random.single.cells.10.nodes.stats, 'single', 'random_10')
compare.performance.plot.2d(experiments.random.single.cells.15.nodes.stats, 'single', 'random_15')
compare.performance.plot.2d(experiments.random.single.cells.20.nodes.stats, 'single', 'random_20')

compare.performance.plot.2d(experiments.random.multiple.biopses.5.nodes.stats, 'multiple', 'random_5')
compare.performance.plot.2d(experiments.random.multiple.biopses.10.nodes.stats, 'multiple', 'random_10')
compare.performance.plot.2d(experiments.random.multiple.biopses.15.nodes.stats, 'multiple', 'random_15')
compare.performance.plot.2d(experiments.random.multiple.biopses.20.nodes.stats, 'multiple', 'random_20')

### end of file -- statistics.random.main.R

