##############################################################################
###
### MST
###
### Test mltree
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

exp = experiments.random.single.cells.5.nodes.scite[[6, 81]][[1]]

data = exp$dataset

dataset = import.genotypes(data)
recon = tronco.mst.edmonds(dataset)
pf = recon$adj.matrix.prima.facie


### end of file -- test.mltree.R
