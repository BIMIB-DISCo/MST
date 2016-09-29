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

load('RData/experiments.random.single.cells.5.nodes.scite.RData')
load('RData/experiments.random.single.cells.10.nodes.scite.RData')
load('RData/experiments.random.single.cells.15.nodes.scite.RData')
load('RData/experiments.random.single.cells.20.nodes.scite.RData')

experiments.random.single.cells.5.nodes.scite.stats = get.stats(experiments.random.single.cells.5.nodes.scite)
save(experiments.random.single.cells.5.nodes.scite.stats, file="RData/experiments.random.single.cells.5.nodes.scite.stats.RData")
giulio.plot(experiments.random.single.cells.5.nodes.scite.stats, 'single', 'random_5')

experiments.random.single.cells.10.nodes.scite.stats = get.stats(experiments.random.single.cells.10.nodes.scite)
save(experiments.random.single.cells.10.nodes.scite.stats, file="RData/experiments.random.single.cells.10.nodes.scite.stats.RData")
giulio.plot(experiments.random.single.cells.10.nodes.scite.stats, 'single', 'random_10')

experiments.random.single.cells.15.nodes.scite.stats = get.stats(experiments.random.single.cells.15.nodes.scite)
save(experiments.random.single.cells.15.nodes.scite.stats, file="RData/experiments.random.single.cells.15.nodes.scite.stats.RData")
giulio.plot(experiments.random.single.cells.15.nodes.scite.stats, 'single', 'random_15')

experiments.random.single.cells.20.nodes.scite.stats = get.stats(experiments.random.single.cells.20.nodes.scite)
save(experiments.random.single.cells.20.nodes.scite.stats, file="RData/experiments.random.single.cells.20.nodes.scite.stats.RData")
giulio.plot(experiments.random.single.cells.20.nodes.scite.stats, 'single', 'random_20')

