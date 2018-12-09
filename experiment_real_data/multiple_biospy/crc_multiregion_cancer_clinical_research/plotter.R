##############################################################################
###
### MST
###
### Plotter
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


library(TRONCO)

f = "CAPRESE-clean.Rdata"
load(f) 
tronco.plot(m,
            scale.nodes = .6,
            confidence = c('tp', 'pr', 'hg'),
            legend.cex = .5,
            label.edge.size = 9,
            legend.pos = 'top',
            disconnected = T)
dev.copy2pdf(file = paste(f,'.pdf', sep = ''))

### end of file -- plotter.R
