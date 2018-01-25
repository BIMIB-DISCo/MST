

source('../generate.sifit.input.R')

load('RData/dataset.multiple.biopses.random.columns.low.RData')
load('RData/dataset.multiple.biopses.random.columns.medium.RData')
load('RData/dataset.multiple.biopses.random.columns.high.RData')


cat('scite medium\n')
create.sifit.input(dataset.multiple.biopses.random.columns.medium, 'multiple', 'medium', scite.sd, numsample = 10, pass.error.rates = FALSE)