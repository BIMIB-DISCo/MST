library(TRONCO)
library(ape)

load('data.RData')
complete = data
geno = as.genotypes(complete)
geno = keysToNames(complete, geno)
geno = geno[-6,]
rownames = paste0('s', 1:nrow(geno))
colnames = colnames(geno)
colnames = gsub(' ', '_', colnames)
colnames(geno) = colnames
stree = nj(dist.gene(geno))
plot(stree)
dev.copy2pdf(file = 'crc.pdf')

load('data.merged.RData')
complete = data.merged
geno = as.genotypes(complete)
geno = keysToNames(complete, geno)
geno = geno[-6,]
rownames = paste0('s', 1:nrow(geno))
colnames = colnames(geno)
colnames = gsub(' ', '_', colnames)
colnames(geno) = colnames
stree = nj(dist.gene(geno))
plot(stree)
dev.copy2pdf(file = 'crc_merged.pdf')
