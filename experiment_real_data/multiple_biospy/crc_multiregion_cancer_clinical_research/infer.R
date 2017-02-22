install.packages('/Volumes/DATA/Work/Software/Github/TRONCO_2.4.2.tar.gz')

library(RColorBrewer)
library(TRONCO)

SKIPDATA = TRUE
DOPLOT = FALSE
NBOOT = 1000

file = "./data/Workbook1.XLSX"

if (!SKIPDATA) {

    install.packages("xlsx", dependencies = T, repo = "http://cran.us.r-project.org")
    install.packages("rJava", dependencies = T, repo = "http://cran.us.r-project.org")
    library(xlsx)

    y = read.xlsx(file, sheetIndex = 1)[1:208, ] #All events
    # y = read.xlsx(file, sheetIndex = 3)[1:18, ] # clean
    types = as.vector(unique(y$Mutation))
        
    to.TRONCO = function(t)
    {
        d = y[y$Mutation == t, ]
        d = d[, c('Gene', 'P1', 'P2', 'P3', 'P4', 'M1', 'M2', 'M3')]
        rownames(d) = d$Gene
        d = d[, c('P1', 'P2', 'P3', 'P4', 'M1', 'M2', 'M3')]
        d = t(d)
        return(import.genotypes(d, t))
    }    

    data = NULL
    for(i in 1:length(types))
    {
        tmp = to.TRONCO(types[i])
        
        data$genotypes = cbind(data$genotypes, tmp$genotypes)
        data$annotations = rbind(data$annotations, tmp$annotations)
        colnames(data$genotypes) = paste('G', 1:ncol(data$genotypes), sep='')
        rownames(data$annotations) = paste('G', 1:nrow(data$annotations), sep='')
        data$types = rbind(data$types, tmp$types)
    }
    
    
    data$types[, 1] = colorRampPalette(brewer.pal(8, "Accent"))(nrow(data$types))
    data = enforce.numeric(data)
    save(data, file = "MAF.Rdata")
}
if (SKIPDATA) { load("MAF.Rdata") }


if(DOPLOT)
{
    dev.new(height=22)
    source('/Volumes/DATA/Work/Software/Github/TRONCO/R/visualization.R', chdir = F)
    oncoprint(data, sample.id = T,  genes.cluster = T,  font.column = 9, cellwidth = 10, font.row = 8, annotate.consolidate.events = F)
    dev.copy2pdf(file = 'data-pat1.pdf')
}


data$genotypes = rbind(data$genotypes, rep(0, ncol(data$genotypes)))
consolidate.data(data)

# CHOWLIU si impalla
# ALGO = c("EDMONDS", "GABOW", "PRIM", "CAPRI", "CAPRESE")
ALGO = c("EDMONDS", "PRIM", "CAPRI",  "CHOWLIU", "CAPRESE")


infer = function(x, y, ann, ...) {
    if (y == "PRIM") 
        m = tronco.mst.prim(x, nboot = NBOOT)
    if (y == "EDMONDS") 
        m = tronco.mst.edmonds(x, nboot = NBOOT)
    if (y == "CHOWLIU") 
        m = tronco.mst.chowliu(x, nboot = NBOOT)
    if (y == "GABOW") 
        m = tronco.mst.gabow(x, nboot = NBOOT)
    if (y == "CAPRESE") 
        m = tronco.caprese(x)
    if (y == "CAPRI") 
        m = tronco.capri(x, nboot = NBOOT)

    if(DOPLOT)
    {
        tronco.plot(m, ...)
        dev.copy2pdf(file = paste("model-", ann, "-", 
            y, ".pdf", sep = ""))
    }
    
    save(m, file=paste(y, '-', ann, '.Rdata', sep =''))
}

sapply(ALGO, FUN = infer, x = data, ann = "clean")
