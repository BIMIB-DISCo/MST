# example of usage of the package phangorn to reconstruct phylogenetic trees
library(phangorn)

# create sample data starting from the ones we use in TRONCO
data = array("missing",c(4,6))
data[1,3]  = "observed"
data[1,6]  = "observed"
data[3,5]  = "observed"
data[4,5]  = "observed"
data = t(data)
rownames(data) = as.character(1:nrow(data))

# create the input data for the phylogeny algorithms
phylo.data = phyDat(data,type="USER",levels=c("missing","observed"))

# infer the best tree given the hamming distance as metric
dm  = dist.hamming(phylo.data)
hamming.distance.phylo = upgma(dm)

# infer the maximum parsimony phylogenetic tree starting from the hamming distance tree
max.parsimony.phylo  <- optim.parsimony(hamming.distance.phylo, phylo.data)
