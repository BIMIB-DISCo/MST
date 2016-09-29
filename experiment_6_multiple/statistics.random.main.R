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

load('RData/experiments.random.multiple.biopses.forest.nodes.scite.RData')

experiments.random.multiple.biopses.forest.nodes.scite.stats = get.stats(experiments.random.multiple.biopses.forest.nodes.scite)
save(experiments.random.multiple.biopses.forest.nodes.scite.stats, file="RData/experiments.random.multiple.biopses.forest.nodes.scite.stats.RData")
giulio.plot(experiments.random.multiple.biopses.forest.nodes.scite.stats, 'multiple', 'random_forest')

