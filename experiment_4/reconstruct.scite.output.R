##############################################################################
###
### MST
###
### Reconstruct SCITE Output
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

## After recon

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/result.single.cells.mini.01.RData')
load('RData/result.single.cells.mini.02.RData')
load('RData/result.single.cells.mini.03.RData')
load('RData/result.single.cells.mini.04.RData')
load('RData/result.single.cells.mini.05.RData')

library(igraph)
library(sna)
library(Rgraphviz)


## Merge tronco results with SCITE

experiments.single.cells.mini.01.scite =
    import.scite.output(result.single.cells.mini.01, 'single', 'mini_01')
save(experiments.single.cells.mini.01.scite,
     file = 'RData/experiments.single.cells.mini.01.scite.RData')

experiments.single.cells.mini.02.scite =
    import.scite.output(result.single.cells.mini.02, 'single', 'mini_02')
save(experiments.single.cells.mini.02.scite,
     file = 'RData/experiments.single.cells.mini.02.scite.RData')

experiments.single.cells.mini.03.scite =
    import.scite.output(result.single.cells.mini.03, 'single', 'mini_03')
save(experiments.single.cells.mini.03.scite,
     file = 'RData/experiments.single.cells.mini.03.scite.RData')

experiments.single.cells.mini.04.scite =
    import.scite.output(result.single.cells.mini.04, 'single', 'mini_04')
save(experiments.single.cells.mini.04.scite,
     file = 'RData/experiments.single.cells.mini.04.scite.RData')

experiments.single.cells.mini.05.scite =
    import.scite.output(result.single.cells.mini.05, 'single', 'mini_05')
save(experiments.single.cells.mini.05.scite,
     file = 'RData/experiments.single.cells.mini.05.scite.RData')

### end of file -- reconstruct.scite.output.R
