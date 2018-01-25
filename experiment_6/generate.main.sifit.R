

source('../generate.sifit.input.R')

load('RData/dataset.random.single.cells.forest.nodes.RData')


cat('sifit random\n')
create.sifit.input(dataset.random.single.cells.forest.nodes, 'single', 'random_forest', scite.sd, numsample = 50)