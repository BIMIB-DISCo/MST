##############################################################################
###
### MST
###
### Analysis
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


## Load the required R packages

library(TRONCO)

## set the seed
set.seed(483775)

epos = 1.24e-6
eneg = 0.0972

## load the data
load("dataset_navin_TNBC.RData")
dataset = import.genotypes(dataset_navin_TNBC)

## infer the progressions
capri = tronco.capri(dataset, epos = epos, eneg = eneg)
caprese = tronco.caprese(dataset, epos = epos, eneg = eneg)
edmonds = tronco.bootstrap(tronco.edmonds(dataset, epos = epos, eneg = eneg))
gabow = tronco.gabow(dataset, epos = epos, eneg = eneg)
chowliu = tronco.chowliu(dataset, epos = epos, eneg = eneg)
prim = tronco.prim(dataset, epos = epos, eneg = eneg)

## save the results
save(capri, file = "capri.RData")
save(caprese, file = "caprese.RData")
save(edmonds, file = "edmonds.RData")
save(gabow, file = "gabow.RData")
save(chowliu, file = "chowliu.RData")
save(prim, file = "prim.RData")

plot_graph <- function(model, title, file) {
    model = change.color(model, type = 'variant', new.color = 'lightblue3')
    oncoprint(model,
              file = paste0('plot/', file, '_data.pdf'),
              title = title)
    tronco.plot(model,
                edge.cex = 1.5,          # scale edge size
                legend.cex = .5,         # scale legend size
                scale.nodes = .6,        # scale node size
                confidence = c('tp', 'pr', 'hg'), # display p-values for these statistics
                disconnected = FALSE,        # do not display nodes without incoming/outgoing edges
                height.logic = .3,       # scale logical connectives
                file = paste0('plot/', file, '.pdf'),
                title = title)
}

plot_graph(capri, title = "CAPRI", file = "capri")
plot_graph(caprese, title = "CAPRESE", file = "caprese")
plot_graph(edmonds, title = "EDMONDS", file = "edmonds")
plot_graph(gabow, title = "GABOW", file = "gabow")
plot_graph(chowliu, title = "CHOW LIU", file = "chowliu")
plot_graph(prim, title = "PRIM", file = "prim")

### end of file -- analysis.R
