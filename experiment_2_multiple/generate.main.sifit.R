

source('../generate.sifit.input.R')

load('RData/dataset.multiple.biopses.low.RData')
load('RData/dataset.multiple.biopses.medium.RData')
load('RData/dataset.multiple.biopses.high.RData')

cat('scite medium\n')
create.sifit.input(dataset.multiple.biopses.medium, 'multiple', 'medium', scite.sd, numsample = 10, pass.error.rates = FALSE)