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


## Set the seed

SEED = 483775

epos = 1.24e-6
eneg = 0.0972


## Load the data

dataset = import.genotypes(dataset_navin_TNBC)


## Infer the progressions

edmonds =
    tronco.bootstrap(tronco.edmonds(dataset,
                                    epos = epos,
                                    eneg = eneg,
                                    boot.seed = SEED))


## Save the results

save(edmonds, file = "edmonds.RData")


plot_graph <- function(model, title, file) {
    model = change.color(model, type = 'variant', new.color = 'lightblue3')
    oncoprint(model, file = paste0(file, '_data.pdf'), title = title)
    tronco.plot(model,
                edge.cex = 1.5,          # scale edge size
                legend.cex = .5,         # scale legend size
                scale.nodes = .6,        # scale node size
                confidence = c('tp', 'pr', 'hg'), # display p - values for these statistics
                disconnected = FALSE,        # do not display nodes without incoming/outgoing edges
                height.logic = .3,       # scale logical connectives
                file = paste0(file, '.pdf'),
                title = title)
}

plot_graph(edmonds, title = "EDMONDS", file = "edmonds")

### end of file -- analysis.R
