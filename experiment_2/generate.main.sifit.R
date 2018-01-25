

source('../generate.sifit.input.R')

load('RData/dataset.single.cells.low.RData')
load('RData/dataset.single.cells.medium.RData')
load('RData/dataset.single.cells.high.RData')

#cat('sifit low\n')
#create.sifit.input(dataset.single.cells.low, 'single', 'low', scite.sd)
cat('sifit medium\n')
create.sifit.input(dataset.single.cells.medium, 'single', 'medium', scite.sd, 50)
#cat('sifit high\n')
#create.sifit.input(dataset.single.cells.high, 'single', 'high', scite.sd)
