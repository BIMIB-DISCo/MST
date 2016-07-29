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

load('RData/result.random.multiple.biopses.bulk.5.nodes.RData')
load('RData/result.random.multiple.biopses.bulk.10.nodes.RData')
load('RData/result.random.multiple.biopses.bulk.15.nodes.RData')
load('RData/result.random.multiple.biopses.bulk.20.nodes.RData')



library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiments.random.multiple.biopses.bulk.5.nodes.scite = import.scite.output(result.random.multiple.biopses.bulk.5.nodes, 'multiple', 'random_5')
save(experiments.random.multiple.biopses.bulk.5.nodes.scite, file = 'RData/experiments.random.multiple.biopses.bulk.5.nodes.scite.RData')

experiments.random.multiple.biopses.bulk.10.nodes.scite = import.scite.output(result.random.multiple.biopses.bulk.10.nodes, 'multiple', 'random_10')
save(experiments.random.multiple.biopses.bulk.10.nodes.scite, file = 'RData/experiments.random.multiple.biopses.bulk.10.nodes.scite.RData')

experiments.random.multiple.biopses.bulk.15.nodes.scite = import.scite.output(result.random.multiple.biopses.bulk.15.nodes, 'multiple', 'random_15')
save(experiments.random.multiple.biopses.bulk.15.nodes.scite, file = 'RData/experiments.random.multiple.biopses.bulk.15.nodes.scite.RData')

experiments.random.multiple.biopses.bulk.20.nodes.scite = import.scite.output(result.random.multiple.biopses.bulk.20.nodes, 'multiple', 'random_20')
save(experiments.random.multiple.biopses.bulk.20.nodes.scite, file = 'RData/experiments.random.multiple.biopses.bulk.20.nodes.scite.RData')
