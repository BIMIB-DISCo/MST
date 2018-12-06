##############################################################################
###
### MST
###
### Generate Main Sifit
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


source('../generate.sifit.input.R')

load('RData/dataset.multiple.biopses.low.RData')
load('RData/dataset.multiple.biopses.medium.RData')
load('RData/dataset.multiple.biopses.high.RData')

cat('scite medium\n')
create.sifit.input(dataset.multiple.biopses.medium,
                   'multiple',
                   'medium',
                   scite.sd,
                   numsample = 10,
                   pass.error.rates = FALSE)

### end of file -- generate.main.sifit.R
