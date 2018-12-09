##############################################################################
###
### MST
###
### Generate pf
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


library(TRONCO)
library(ape)

load('data.RData')
complete = data
geno = as.genotypes(complete)
geno = keysToNames(complete, geno)
geno = geno[-6, ]
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
geno = geno[-6, ]
rownames = paste0('s', 1:nrow(geno))
colnames = colnames(geno)
colnames = gsub(' ', '_', colnames)
colnames(geno) = colnames
stree = nj(dist.gene(geno))
plot(stree)
dev.copy2pdf(file = 'crc_merged.pdf')

### end of file -- generate.pf.R
