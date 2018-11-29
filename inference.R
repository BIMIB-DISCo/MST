##############################################################################
###
### MST
###
### Inference
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


library( devtools )  
install_github( "BIMIB-DISCo/TRONCO", ref = 'development' )


SKIPDATA = TRUE
DOPLOT = TRUE
NBOOT = 100

primary = paste("p3-", 1:3, sep = "")
metastatic = paste("L", 1:2, sep = "")
regions = c(primary, metastatic)

file = "./data/journal.pone.0152673.s002.XLSX"

if (!SKIPDATA) {

    install.packages("xlsx", dependencies = T, repo = "http://cran.us.r-project.org")
    install.packages("rJava", dependencies = T, repo = "http://cran.us.r-project.org")
    library(xlsx)

    to.MAF = function(x) {
        y = read.xlsx(file, sheetIndex = x)[, c("Gene.knownGene", 
                                                "ExonicFunc.knownGene")]
        y$Tumor_Sample_Barcode = rep(x, nrow(y))
        colnames(y) = c("Hugo_Symbol", "Variant_Classification", 
			"Tumor_Sample_Barcode")
        return(y)
    }

    p3.1 = to.MAF(regions[1])
    p3.2 = to.MAF(regions[2])
    p3.3 = to.MAF(regions[3])
    L1 = to.MAF(regions[4])
    L2 = to.MAF(regions[5])

    MAF = rbind(p3.1, p3.2, p3.3, L1, L2)
    save(MAF, file = "MAF.Rdata")
}


if (SKIPDATA) { load("MAF.Rdata") }

## CHOWLIU si impalla (agosto 2016)
## ALGO = c("EDMONDS", "GABOW", "PRIM", "CAPRI", "CAPRESE")

ALGO = c("EDMONDS", "PRIM", "CAPRI", "CAPRESE")

infer <- function(x, y, ann, ...) {
    if (y == "PRIM") 
        m = tronco.prim(x, nboot = NBOOT)
    if (y == "EDMONDS") 
        m = tronco.edmonds(x, nboot = NBOOT)
    if (y == "CHOWLIU") 
        m = tronco.chowliu(x, nboot = NBOOT)
    if (y == "GABOW") 
        m = tronco.gabow(x, nboot = NBOOT)
    if (y == "CAPRESE") 
        m = tronco.caprese(x)
    if (y == "CAPRI") 
        m = tronco.capri(x, nboot = NBOOT)

    m = tronco.bootstrap(m, nboot = NBOOT)
    
    if (DOPLOT)
    {
        tronco.plot(m, ...)
        dev.copy2pdf(file = paste("model-", ann, "-", 
                                  y, ".pdf", sep = ""))
    }
    
    save(m, file=paste(y, '-', ann, '.Rdata', sep =''))
}


library(TRONCO)


## Model with all events (undistinguishable + 1 fake sample for p(x) < 1)
data = import.MAF(MAF, merge.mutation.types = FALSE)
data$genotypes = rbind(data$genotypes, rep(0, ncol(data$genotypes)))

tofix = consolidate.data(data, print = T)

i f(DOPLOT) {
    oncoprint(data, sample.id = T, font.column = 9, cellwidth = 10, genes.cluster = T)
    dev.copy2pdf(file = "oncoprint.pdf")
}

sapply(ALGO, FUN = infer, x = data, ann = "")

## Model with all clean events (1 fake sample for p(x) < 1)
data = import.MAF(MAF, merge.mutation.types = FALSE)
err = consolidate.data(data)

## Delete all 1s
for (i in 1:length(err$ones)) {
    errs = err$ones[[i]]
    data = delete.event(data, errs["event"], errs["type"])
}

## Add a fake event with all 1s, plus a fake sample
data$genotypes = cbind(data$genotypes, rep(1, nrow(data$genotypes))) # event
colnames(data$genotypes)[ncol(data$genotypes)] = "M"

data$annotations = rbind(data$annotations,
                         c("indistinguishable", "clonal"))
rownames(data$annotations)[nrow(data$annotations)] = "M"

data$types = rbind(data$types, c("brown3"))
rownames(data$types)[nrow(data$types)] = "indistinguishable"

data

data$genotypes = rbind(data$genotypes,
                       rep(0, ncol(data$genotypes))) # fake sample


## Merge indistinguishable G37 and G38 -- change name and type to G39
tofix = consolidate.data(data, print = T)
data$genotypes = data$genotypes[, which(colnames(data$genotypes) !=
                                        "G37")]
data$genotypes = data$genotypes[, which(colnames(data$genotypes) != 
                                        "G38")]
data$annotations = data$annotations[which(rownames(data$annotations) != 
                                          "G37"), ]
data$annotations = data$annotations[which(rownames(data$annotations) != 
                                          "G38"), ]
data$annotations["G39", "event"] = "merged"
data$annotations["G39", "type"] = "indistinguishable"
data$types = data$types[which(rownames(data$types) != "NA"),
                       ,
                        drop = FALSE]

if (DOPLOT) {
    oncoprint(data, sample.id = T, font.column = 9, cellwidth = 10)
    dev.copy2pdf(file = "oncoprint-ind.pdf")
}

## ALGO = c("EDMONDS", "GABOW", "PRIM", "CHOWLIU", "CAPRI", "CAPRESE")
ALGO = c("EDMONDS", "PRIM", "GABOW", "CHOWLIU", "CAPRI", "CAPRESE")


## Use these data to execute SCITE

sapply(ALGO,
       FUN = infer,
       x = data,
       ann = "ind",
       scale.nodes = 0.2,
       confidence = c('tp', 'pr', 'hg', 'npb'))

## CAPRI with homologous with all clean events (1 fake sample for p(x) < 1)
data = import.MAF(MAF, merge.mutation.types = FALSE)
err = consolidate.data(data)

data = hypothesis.add.homologous(data)
sapply("CAPRI", FUN = infer, x = data, ann = "hypo", 
	scale.nodes = 0.2)


### end of file -- inference.R
