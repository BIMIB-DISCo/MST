### after recon

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/experiment.missing.data.RData')

library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiment.missing.data.scite = import.scite.output(experiment.missing.data, 'single', 'missing')
save(experiment.missing.data.scite, file = 'RData/experiment.missing.data.scite.RData')
