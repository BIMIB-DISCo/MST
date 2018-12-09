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

epos_level = 0.05
eneg_level = 0.05

load('dataset_navin_TNBC.RData')
data = dataset_navin_TNBC
onconem.genotypes = t(data)

oNEM = oncoNEM$new(Data = onconem.genotypes,
                   FPR = epos_level,
                   FNR = eneg_level)

oNEM$search(delta = 200)

plotTree(tree = oNEM$best$tree,
         clones = NULL,
         vertex.size = 25)
dev.copy2pdf(file = 'onconem.best.pdf')

oNEM.expanded = expandOncoNEM(oNEM,
                              epsilon = 10,
                              delta = 200,
                              checkMax = 10000,
                              app = TRUE)

plotTree(tree = oNEM.expanded$best$tree,
         clones = NULL,
         vertex.size = 25)
dev.copy2pdf(file = 'onconem.best.expand.pdf')

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

dev.copy2pdf(file = 'onconem.best.clustered.pdf')

### end of file -- generate.onconem.R
