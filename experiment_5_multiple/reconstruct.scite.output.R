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

load('RData/result.multiple.biopses.low.RData')
load('RData/result.multiple.biopses.medium.RData')
load('RData/result.multiple.biopses.high.RData')

library(igraph)
library(sna)
library(Rgraphviz)


## Merge tronco results with scite

experiments.multiple.biopses.low.scite =
    import.scite.output(result.multiple.biopses.low, 'multiple', 'low')
save(experiments.multiple.biopses.low.scite,
     file = 'RData/experiments.multiple.biopses.low.scite.RData')

experiments.multiple.biopses.medium.scite =
    import.scite.output(result.multiple.biopses.medium, 'multiple', 'medium')
save(experiments.multiple.biopses.medium.scite,
     file = 'RData/experiments.multiple.biopses.medium.scite.RData')

experiments.multiple.biopses.high.scite =
    import.scite.output(result.multiple.biopses.high, 'multiple', 'high')
save(experiments.multiple.biopses.high.scite,
     file = 'RData/experiments.multiple.biopses.high.scite.RData')

### end of file -- reconstruct.scite.output.R
