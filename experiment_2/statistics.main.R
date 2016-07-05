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

load('RData/experiments.single.cells.low.scite.RData')
load('RData/experiments.single.cells.medium.scite.RData')
load('RData/experiments.single.cells.high.scite.RData')

experiments.single.cells.low.scite.stats = get.stats(experiments.single.cells.low.scite)
save(experiments.single.cells.low.scite.stats, file="RData/experiments.single.cells.low.scite.stats.RData")

experiments.single.cells.medium.scite.stats = get.stats(experiments.single.cells.medium.scite)
save(experiments.single.cells.medium.scite.stats, file="RData/experiments.single.cells.medium.scite.stats.RData")

experiments.single.cells.high.scite.stats = get.stats(experiments.single.cells.high.scite)
save(experiments.single.cells.high.scite.stats, file="RData/experiments.single.cells.high.scite.stats.RData")

