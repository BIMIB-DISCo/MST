##################################################################################
#                                                                                #
# MST                                                                            #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Giulio Caravagna, Luca De Sano, Daniele Ramazzotti         #
# email: tronco@disco.unimib.it                                                  #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################

# source the needed script

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/result.single.cells.clean.RData')
load('RData/result.single.cells.convergent.RData')
load('RData/result.single.cells.random.columns.RData')
load('RData/result.single.cells.random.forest.RData')

library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiments.single.cells.clean.scite = import.scite.output(result.single.cells.clean, 'single', 'clean', c(75))
save(experiments.single.cells.clean.scite, file = 'RData/experiments.single.cells.clean.scite.RData')

experiments.single.cells.convergent.scite = import.scite.output(result.single.cells.convergent, 'single', 'convergent', c(75))
save(experiments.single.cells.convergent.scite, file = 'RData/experiments.single.cells.convergent.scite.RData')

experiments.single.cells.random.columns.scite = import.scite.output(result.single.cells.random.columns, 'single', 'random_columns', c(75))
save(experiments.single.cells.random.columns.scite, file = 'RData/experiments.single.cells.random.columns.scite.RData')

experiments.single.cells.random.forest.scite = import.scite.output(result.single.cells.random.forest, 'single', 'random_forest_fixed', c(75))
save(experiments.single.cells.random.forest.scite, file = 'RData/experiments.single.cells.random.forest.scite.RData')
