##############################################################################
###
### MST
###
### Main
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


if (!require(TRONCO)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("TRONCO")
    library(TRONCO)
}


load('data.RData')

NBOOT = 100
ALGO = c("EDMONDS", "PRIM")
DOPLOT = TRUE

infer <- function(x, y, ann, ...) {
    if (y == "PRIM") 
        m = tronco.prim(x, nboot = NBOOT, boot.seed = 674523)
    if (y == "EDMONDS") 
        m = tronco.edmonds(x, nboot = NBOOT, boot.seed = 10000)
    if (y == "CHOWLIU") 
        m = tronco.chowliu(x, nboot = NBOOT, boot.seed = 234354)
    if (y == "GABOW") 
        m = tronco.gabow(x, nboot = NBOOT, boot.seed = 112243)
    if (y == "CAPRESE") 
        m = tronco.caprese(x)
    if (y == "CAPRI") 
        m = tronco.capri(x, nboot = NBOOT, boot.seed = 20000)

    m = tronco.bootstrap(m, nboot = NBOOT)
    
    if (DOPLOT) {
        tronco.plot(m, ...)
        dev.copy2pdf(file = paste("model-", ann, "-", y, ".pdf",
                                  sep = ""))
    }
    save(m, file=paste(y, '-', ann, '.Rdata', sep =''))
}


if (DOPLOT) {
    oncoprint(data, sample.id = TRUE, font.column = 9, cellwidth = 10)
    dev.copy2pdf(file = "oncoprint-ind.pdf")
}


sapply(ALGO,
       FUN = infer,
       x = data,
       ann = "ind",
       scale.nodes = 0.2,
       confidence = c('tp', 'pr', 'hg', 'npb'))

### end of file -- main.R
