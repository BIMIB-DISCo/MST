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


## Statistics

source('../statistics.plot.R')
source('../statistics.compute.R')
source('../giulio.plot.R')

load('RData/experiments.single.cells.mini.01.scite.RData')
load('RData/experiments.single.cells.mini.02.scite.RData')
load('RData/experiments.single.cells.mini.03.scite.RData')
load('RData/experiments.single.cells.mini.04.scite.RData')
load('RData/experiments.single.cells.mini.05.scite.RData')

experiments.single.cells.mini.01.scite.stats =
    get.stats(experiments.single.cells.mini.01.scite)
save(experiments.single.cells.mini.01.scite.stats,
     file = "RData/experiments.single.cells.mini.01.scite.stats.RData")
giulio.plot(experiments.single.cells.mini.01.scite.stats, 'single', 'mini_01')

experiments.single.cells.mini.02.scite.stats =
    get.stats(experiments.single.cells.mini.02.scite)
save(experiments.single.cells.mini.02.scite.stats,
     file = "RData/experiments.single.cells.mini.02.scite.stats.RData")
giulio.plot(experiments.single.cells.mini.02.scite.stats, 'single', 'mini_02')

experiments.single.cells.mini.03.scite.stats =
    get.stats(experiments.single.cells.mini.03.scite)
save(experiments.single.cells.mini.03.scite.stats,
     file = "RData/experiments.single.cells.mini.03.scite.stats.RData")
giulio.plot(experiments.single.cells.mini.03.scite.stats, 'single', 'mini_03')

experiments.single.cells.mini.04.scite.stats =
    get.stats(experiments.single.cells.mini.04.scite)
save(experiments.single.cells.mini.04.scite.stats,
     file = "RData/experiments.single.cells.mini.04.scite.stats.RData")
giulio.plot(experiments.single.cells.mini.04.scite.stats, 'single', 'mini_04')

experiments.single.cells.mini.05.scite.stats =
    get.stats(experiments.single.cells.mini.05.scite)
save(experiments.single.cells.mini.05.scite.stats,
     file = "RData/experiments.single.cells.mini.05.scite.stats.RData")
giulio.plot(experiments.single.cells.mini.05.scite.stats, 'single', 'mini_05')

### end of file -- statistics.main.R
