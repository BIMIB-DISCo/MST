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

load('RData/result.multiple.biopses.clean.RData')
load('RData/result.multiple.biopses.convergent.RData')
load('RData/result.multiple.biopses.random.columns.RData')
load('RData/result.multiple.biopses.random.forest.RData')

library(igraph)
library(sna)
library(Rgraphviz)


### Merge tronco results with scite

experiments.multiple.biopses.clean.scite =
    import.scite.output(result.multiple.biopses.clean,
                        'multiple',
                        'clean',
                        c(20))
save(experiments.multiple.biopses.clean.scite,
     file = 'RData/experiments.multiple.biopses.clean.scite.RData')


experiments.multiple.biopses.convergent.scite =
    import.scite.output(result.multiple.biopses.convergent,
                        'multiple',
                        'convergent',
                        c(20))
save(experiments.multiple.biopses.convergent.scite,
     file = 'RData/experiments.multiple.biopses.convergent.scite.RData')


experiments.multiple.biopses.random.columns.scite =
    import.scite.output(result.multiple.biopses.random.columns,
                        'multiple',
                        'random_columns',
                        c(20))
save(experiments.multiple.biopses.random.columns.scite,
     file = 'RData/experiments.multiple.biopses.random.columns.scite.RData')


experiments.multiple.biopses.random.forest.scite =
    import.scite.output(result.multiple.biopses.random.forest,
                        'multiple',
                        'random_forest_fixed',
                        c(20))
save(experiments.multiple.biopses.random.forest.scite,
     file = 'RData/experiments.multiple.biopses.random.forest.scite.RData')

### end of file -- reconstruct.scite.output.R
