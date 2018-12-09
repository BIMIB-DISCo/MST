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

## Source the needed script

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/result.single.cells.low.RData')
load('RData/result.single.cells.medium.RData')
load('RData/result.single.cells.high.RData')

library(igraph)
library(sna)
library(Rgraphviz)


## Merge tronco results with scite

experiments.single.cells.low.scite =
    import.scite.output(result.single.cells.low, 'single', 'low')
save(experiments.single.cells.low.scite,
     file = 'RData/experiments.single.cells.low.scite.RData')

experiments.single.cells.medium.scite =
    import.scite.output(result.single.cells.medium, 'single', 'medium')
save(experiments.single.cells.medium.scite,
     file = 'RData/experiments.single.cells.medium.scite.RData')

experiments.single.cells.high.scite =
    import.scite.output(result.single.cells.high, 'single', 'high')
save(experiments.single.cells.high.scite,
     file = 'RData/experiments.single.cells.high.scite.RData')

### end of file -- reconstruct.scite.output.R
