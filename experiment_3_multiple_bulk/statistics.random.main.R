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

load('RData/experiments.random.multiple.biopses.bulk.5.nodes.scite.RData')
load('RData/experiments.random.multiple.biopses.bulk.10.nodes.scite.RData')
load('RData/experiments.random.multiple.biopses.bulk.15.nodes.scite.RData')
load('RData/experiments.random.multiple.biopses.bulk.20.nodes.scite.RData')

experiments.random.multiple.biopses.bulk.5.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.bulk.5.nodes.scite)
save(experiments.random.multiple.biopses.bulk.5.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.bulk.5.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.bulk.5.nodes.scite.stats, 'multiple', 'random_5')

experiments.random.multiple.biopses.bulk.10.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.bulk.10.nodes.scite)
save(experiments.random.multiple.biopses.bulk.10.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.bulk.10.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.bulk.10.nodes.scite.stats, 'multiple', 'random_10')

experiments.random.multiple.biopses.bulk.15.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.bulk.15.nodes.scite)
save(experiments.random.multiple.biopses.bulk.15.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.bulk.15.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.bulk.15.nodes.scite.stats, 'multiple', 'random_15')

experiments.random.multiple.biopses.bulk.20.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.bulk.20.nodes.scite)
save(experiments.random.multiple.biopses.bulk.20.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.bulk.20.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.bulk.20.nodes.scite.stats, 'multiple', 'random_20')
