### after recon

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/result.single.cells.mini.01.RData')
load('RData/result.single.cells.mini.02.RData')
load('RData/result.single.cells.mini.03.RData')
load('RData/result.single.cells.mini.04.RData')
load('RData/result.single.cells.mini.05.RData')

library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiments.single.cells.mini.01.scite = import.scite.output(result.single.cells.mini.01, 'single', 'mini_01')
save(experiments.single.cells.mini.01.scite, file = 'RData/experiments.single.cells.mini.01.scite.RData')

experiments.single.cells.mini.02.scite = import.scite.output(result.single.cells.mini.02, 'single', 'mini_02')
save(experiments.single.cells.mini.02.scite, file = 'RData/experiments.single.cells.mini.02.scite.RData')

experiments.single.cells.mini.03.scite = import.scite.output(result.single.cells.mini.03, 'single', 'mini_03')
save(experiments.single.cells.mini.03.scite, file = 'RData/experiments.single.cells.mini.03.scite.RData')

experiments.single.cells.mini.04.scite = import.scite.output(result.single.cells.mini.04, 'single', 'mini_04')
save(experiments.single.cells.mini.04.scite, file = 'RData/experiments.single.cells.mini.04.scite.RData')

experiments.single.cells.mini.05.scite = import.scite.output(result.single.cells.mini.05, 'single', 'mini_05')
save(experiments.single.cells.mini.05.scite, file = 'RData/experiments.single.cells.mini.05.scite.RData')
