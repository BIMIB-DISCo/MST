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

source('../statistics.plot.R')
source('../statistics.compute.R')
source('../giulio.plot.R')

load('RData/experiments.single.cells.clean.scite.RData')
load('RData/experiments.single.cells.convergent.scite.RData')
load('RData/experiments.single.cells.random.columns.scite.RData')
load('RData/experiments.single.cells.random.forest.scite.RData')


experiments.single.cells.clean.scite.stats = get.stats(experiments.single.cells.clean.scite)
save(experiments.single.cells.clean.scite.stats, file="RData/experiments.single.cells.clean.scite.stats.RData")
giulio.plot(experiments.single.cells.clean.scite.stats, 'single', 'clean', c(75))

experiments.single.cells.convergent.scite.stats = get.stats(experiments.single.cells.convergent.scite)
save(experiments.single.cells.convergent.scite.stats, file="RData/experiments.single.cells.convergent.scite.stats.RData")
giulio.plot(experiments.single.cells.convergent.scite.stats, 'single', 'convergent', c(75))

experiments.single.cells.random.columns.scite.stats = get.stats(experiments.single.cells.random.columns.scite)
save(experiments.single.cells.random.columns.scite.stats, file="RData/experiments.single.cells.random.columns.scite.stats.RData")
giulio.plot(experiments.single.cells.random.columns.scite.stats, 'single', 'random_columns', c(75))

experiments.single.cells.random.forest.scite.stats = get.stats(experiments.single.cells.random.forest.scite)
save(experiments.single.cells.random.forest.scite.stats, file="RData/experiments.single.cells.random.forest.scite.stats.RData")
giulio.plot(experiments.single.cells.random.forest.scite.stats, 'single', 'random_forest', c(75))

