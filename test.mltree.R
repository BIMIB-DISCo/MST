

exp = experiments.random.single.cells.5.nodes.scite[[6,81]][[1]]

data = exp$dataset

dataset = import.genotypes(data)
recon = tronco.mst.edmonds(dataset)
pf = recon$adj.matrix.prima.facie
