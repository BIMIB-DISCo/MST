library(TRONCO)

load("CAPRESE-ind.Rdata") 
tronco.plot(m, scale.nodes = .6, confidence = c('tp', 'pr', 'hg'), legend.cex = .5, label.edge.size =9, legend.pos = 'top', disconnected = T)


dev.copy2pdf(file = '')

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

library(TRONCO)

## Model with all events (undistinguishable + 1 fake sample for p(x) < 1)
data = import.MAF(MAF, merge.mutation.types = FALSE)
data$genotypes = rbind(data$genotypes, rep(0, ncol(data$genotypes)))

tofix = consolidate.data(data, print = T)
if(DOPLOT)
{
    oncoprint(data, sample.id = T, font.column = 9, cellwidth = 10)
    dev.copy2pdf(file = "oncoprint.pdf")
}

sapply(ALGO, FUN = infer, x = data, ann = "")

## Model with all clean events (1 fake sample for p(x) < 1)
data = import.MAF(MAF, merge.mutation.types = FALSE)
err = consolidate.data(data)

# delete all 1s
for (i in 1:length(err$ones)) {
    errs = err$ones[[i]]
    data = delete.event(data, errs["event"], errs["type"])
}

# add a fake event with all 1s, plus a fake sample
data$genotypes = cbind(data$genotypes, rep(1, nrow(data$genotypes))) # event
colnames(data$genotypes)[ncol(data$genotypes)] = "M"

data$annotations = rbind(data$annotations, c("indistinguishable", 
    "clonal"))
rownames(data$annotations)[nrow(data$annotations)] = "M"

data$types = rbind(data$types, c("brown3"))
rownames(data$types)[nrow(data$types)] = "indistinguishable"

data

data$genotypes = rbind(data$genotypes, rep(0, ncol(data$genotypes))) # fake sample

# merge indistinguishable G37 and G38 -- change name and type to G39
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
data$types = data$types[which(rownames(data$types) != 
    "NA"), , drop = FALSE]

if(DOPLOT) {
oncoprint(data, sample.id = T, font.column = 9, cellwidth = 10)
dev.copy2pdf(file = "oncoprint-ind.pdf")
}

#ALGO = c("EDMONDS", "GABOW", "PRIM", "CHOWLIU", "CAPRI",     "CAPRESE")
ALGO = c("EDMONDS", "PRIM", "CHOWLIU", "CAPRI", "CAPRESE")

sapply(ALGO, FUN = infer, x = data, ann = "ind", scale.nodes = 0.2)

## CAPRI with homologous with all clean events (1 fake sample for p(x) < 1)
data = import.MAF(MAF, merge.mutation.types = FALSE)
err = consolidate.data(data)

data = hypothesis.add.homologous(data)
sapply("CAPRI", FUN = infer, x = data, ann = "hypo", 
    scale.nodes = 0.2)
