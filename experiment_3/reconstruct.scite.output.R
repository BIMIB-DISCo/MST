##############################################################################
###
### MST
###
### Reconstruct Scite Output
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

load('RData/result.random.single.cells.5.nodes.RData')
load('RData/result.random.single.cells.10.nodes.RData')
load('RData/result.random.single.cells.15.nodes.RData')
load('RData/result.random.single.cells.20.nodes.RData')



library(igraph)
library(sna)
library(Rgraphviz)

### merge tronco results with scite
experiments.random.single.cells.5.nodes.scite = import.scite.output(result.random.single.cells.5.nodes, 'single', 'random_5')
save(experiments.random.single.cells.5.nodes.scite, file = 'RData/experiments.random.single.cells.5.nodes.scite.RData')

experiments.random.single.cells.10.nodes.scite = import.scite.output(result.random.single.cells.10.nodes, 'single', 'random_10')
save(experiments.random.single.cells.10.nodes.scite, file = 'RData/experiments.random.single.cells.10.nodes.scite.RData')

experiments.random.single.cells.15.nodes.scite = import.scite.output(result.random.single.cells.15.nodes, 'single', 'random_15')
save(experiments.random.single.cells.15.nodes.scite, file = 'RData/experiments.random.single.cells.15.nodes.scite.RData')

experiments.random.single.cells.20.nodes.scite = import.scite.output(result.random.single.cells.20.nodes, 'single', 'random_20')
save(experiments.random.single.cells.20.nodes.scite, file = 'RData/experiments.random.single.cells.20.nodes.scite.RData')
