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

load('RData/result.random.single.cells.forest.nodes.RData')


library(igraph)
library(sna)
library(Rgraphviz)

## Merge tronco results with scite

experiments.random.single.cells.forest.nodes.scite =
    import.scite.output(result.random.single.cells.forest.nodes,
                        'single',
                        'random_forest')
save(experiments.random.single.cells.forest.nodes.scite,
     file = 'RData/experiments.random.single.cells.forest.nodes.scite.RData')

### end of file -- reconstruct.scite.output.R
