library(TRONCO)

f = "CAPRESE-clean.Rdata"
 load(f) 
tronco.plot(m, scale.nodes = .6, confidence = c('tp', 'pr', 'hg'), legend.cex = .5, label.edge.size =9, legend.pos = 'top', disconnected = T)
dev.copy2pdf(file=paste(f,'.pdf', sep = ''))
