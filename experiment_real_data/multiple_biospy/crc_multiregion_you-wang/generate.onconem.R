##############################################################################
###
### MST
###
### Generate OncoNEM
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

library(oncoNEM)
library(TRONCO)
options(scipen = 999)

if (! dir.exists('onconem')) {
    dir.create('onconem')
}

epos_level = 0.05
eneg_level = 0.05

load('data.RData')
complete = data
geno = as.genotypes(complete)
geno = keysToNames(complete, geno)
geno = geno[-6, ]
seed_level = round(runif(1) * 10000, 0)
onconem.genotypes = t(geno)

oNEM = oncoNEM$new(Data = onconem.genotypes,
                   FPR = epos_level,
                   FNR = eneg_level)

oNEM$search(delta = 200)

plotTree(tree = oNEM$best$tree, clones = NULL, vertex.size = 25)
dev.copy2pdf(file = 'onconem/onconem.best.pdf')

oNEM.expanded = expandOncoNEM(oNEM,
                              epsilon = 10,
                              delta = 200,
                              checkMax = 10000,
                              app = TRUE)

plotTree(tree = oNEM.expanded$best$tree,
         clones = NULL,
         vertex.size = 25)
dev.copy2pdf(file = 'onconem/onconem.best.expand.pdf')

oncoTree = clusterOncoNEM(oNEM = oNEM.expanded, epsilon = 10)

post = oncoNEMposteriors(tree = oncoTree$g,
                         clones = oncoTree$clones,
                         Data = oNEM$Data,
                         FPR = oNEM$FPR,
                         FNR = oNEM$FNR)

edgeLengths = colSums(post$p_theta)[-1]

plotTree(tree = oncoTree$g,
         clones = oncoTree$clones,
         e.length = edgeLengths,
         label.length = 4,
         axis = TRUE)

dev.copy2pdf(file = 'onconem/onconem.best.clustered.pdf')


load('data.merged.RData')
ind = data.merged
geno = as.genotypes(ind)
geno = keysToNames(ind, geno)
geno = geno[-6, ]
seed_level = round(runif(1) * 10000, 0)
onconem.genotypes = t(geno)

oNEM = oncoNEM$new(Data = onconem.genotypes,
                   FPR = epos_level,
                   FNR = eneg_level)

oNEM$search(delta = 200)

plotTree(tree = oNEM$best$tree, clones = NULL, vertex.size = 25)
dev.copy2pdf(file = 'onconem/onconem.best.merged.pdf')

oNEM.expanded = expandOncoNEM(oNEM,
                              epsilon = 10,
                              delta = 200,
                              checkMax = 10000,
                              app = TRUE)

plotTree(tree = oNEM.expanded$best$tree,
         clones = NULL,
         vertex.size = 25)
dev.copy2pdf(file = 'onconem/onconem.best.merged.expand.pdf')

oncoTree = clusterOncoNEM(oNEM = oNEM.expanded, epsilon = 10)

post = oncoNEMposteriors(tree = oncoTree$g,
                         clones = oncoTree$clones,
                         Data = oNEM$Data,
                         FPR = oNEM$FPR,
                         FNR = oNEM$FNR)

edgeLengths = colSums(post$p_theta)[-1]

plotTree(tree = oncoTree$g,
         clones = oncoTree$clones,
         e.length = edgeLengths,
         label.length = 4,
         axis = TRUE)

dev.copy2pdf(file = 'onconem/onconem.best.merged.clustered.pdf')

### end of file -- generate.onconem.R
