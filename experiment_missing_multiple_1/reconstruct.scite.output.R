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

## After recon

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/experiment.missing.data.RData')

library(igraph)
library(sna)
library(Rgraphviz)


## Merge TRONCO results with SCITE

experiment.missing.data.scite =
    import.scite.output(experiment.missing.data, 'single', 'missing')
save(experiment.missing.data.scite,
     file = 'RData/experiment.missing.data.scite.RData')

### end of file -- reconstruct.scite.output.R
