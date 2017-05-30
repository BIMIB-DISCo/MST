library(TRONCO)

load('data.RData')

NBOOT = 100
ALGO = c("EDMONDS", "PRIM")
SEED = 33333
DOPLOT = TRUE

infer = function(x, y, ann, ...) {
    if (y == "PRIM") 
        m = tronco.prim(x, nboot = NBOOT, boot.seed = SEED)
    if (y == "EDMONDS") 
        m = tronco.edmonds(x, nboot = NBOOT, boot.seed = SEED)
    if (y == "CHOWLIU") 
        m = tronco.chowliu(x, nboot = NBOOT, boot.seed = SEED)
    if (y == "GABOW") 
        m = tronco.gabow(x, nboot = NBOOT, boot.seed = SEED)
    if (y == "CAPRESE") 
        m = tronco.caprese(x)
    if (y == "CAPRI") 
        m = tronco.capri(x, nboot = NBOOT, boot.seed = SEED)

    m = tronco.bootstrap(m, nboot = NBOOT)
    
    if(DOPLOT)
    {
        tronco.plot(m, ...)
        dev.copy2pdf(file = paste("model-", ann, "-", 
            y, ".pdf", sep = ""))
    }
    
    save(m, file=paste(y, '-', ann, '.Rdata', sep =''))
}

if(DOPLOT) {
	oncoprint(data, sample.id = TRUE, font.column = 9, cellwidth = 10)
	dev.copy2pdf(file = "oncoprint-ind.pdf")
}

sapply(ALGO, FUN = infer, x = data, ann = "ind", scale.nodes = 0.2, confidence = c('tp', 'pr', 'hg', 'npb'))
