##############################################################################
###
### MST
###
### Experiment OncoNEM
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

load('RData/dataset.single.cells.medium.RData')

dataset = dataset.single.cells.medium["75", , drop = FALSE]

for (i in 1:ncol(dataset)) {
    exp = dataset[[1, i]][[2]]
    exp = list("1" = exp)
    dataset[[1, i]] = exp
}

start = proc.time()
total = c()
for (i in 1:100) {
    data = dataset[[1, i]][[1]]

    epos = data$epos
    eneg = data$eneg
    genotype = data$dataset

    a = proc.time()

    tronco = import.genotypes(genotype)
    tronco = tronco.gabow(tronco, epos = epos, eneg = eneg)

    b = proc.time() - a
    time = b[['elapsed']]
    tronco.plot(tronco, title = paste("tronco capri time: ", time))
    dev.copy2pdf(file = paste0("plot/", i, "_gabow.pdf"))
    print(time)
    total = c(total, time)

}


stop = proc.time() - start
print(stop[['elapsed']])
mean(total)

create.scite.input(dataset, "single", "medium", 0)

for (i in 1:15) {
    data = dataset[[1, i]][[1]]

    epos = data$epos
    eneg = data$eneg
    onconem.dataset = t(data$dataset)

    a = proc.time()

    recon.onconem = oncoNEM$new(Data = onconem.dataset,
                                FPR = epos,
                                FNR = eneg)
    recon.onconem$search(delta = 50, verbose = TRUE)
    ## linear.tree.onconem = recon.onconem$best$tree
    score.onconem = recon.onconem$best$llh
    oncoTree = clusterOncoNEM(oNEM = recon.onconem, epsilon = 10)
    post = oncoNEMposteriors(tree = oncoTree$g,
                             clones = oncoTree$clones,
                             Data = recon.onconem$Data,
                             FPR = recon.onconem$FPR,
                             FNR = recon.onconem$FNR)
    edgeLengths = colSums(post$p_theta)[-1]

    dev.new()
    plotTree(tree = oncoTree$g,
             clones = oncoTree$clones,
             e.length = edgeLengths,
             label.length = 4)

    b = proc.time() - a
    time = b[['elapsed']]
    print(time)

    title(paste('\noncoNEM ', round(score.onconem, 3), "time: ", time))
    dev.copy2pdf(file = paste0("plot/", i, '_onconem.pdf'))
    dev.off()
}

### end of file -- experiment.R
