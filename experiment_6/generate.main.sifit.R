##############################################################################
###
### MST
###
### Generate Main SIFIT
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


source('../generate.sifit.input.R')

load('RData/dataset.random.single.cells.forest.nodes.RData')


cat('sifit random\n')
create.sifit.input(dataset.random.single.cells.forest.nodes,
                   'single',
                   'random_forest',
                   scite.sd,
                   numsample = 50)

### end of file -- generate.main.sifit.R
