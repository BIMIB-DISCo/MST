# load library
library(igraph)
library(bnlearn)
library(Rgraphviz)
library(sna)
library(oncoNEM)
library(TRONCO)


names = NULL

load('CAPRESE-.Rdata')
complete = m
geno = as.genotypes(complete)
geno = keysToNames(complete, geno)
names$complete = colnames(geno)

load('CAPRESE-ind.Rdata')
ind = m
geno = as.genotypes(ind)
geno = keysToNames(ind, geno)
names$ind = colnames(geno)



for (name in c('complete', 'ind')) {
    filename = paste0('scite_output/datasets/', name, '_ml0')
    read = readLines(paste0(filename, '.gv'))
    read = gsub(' ', '', read)
    read = gsub(';', '', read)
    write(paste0(read, collapse = '\n'), file = paste0(filename, '.correct.gv'))
    readdot = read.dot(paste0(filename, '.correct.gv'))
    graph = graph.adjacency(readdot)
    scite.tree.raw = get.adjacency(graph, sparse = FALSE)
    scite.tree = matrix(0, ncol = ncol(scite.tree.raw) - 1, nrow = ncol(scite.tree.raw) - 1)
    for (i in 1:nrow(scite.tree)) {
        for(j in 1:ncol(scite.tree)) {
            if (scite.tree.raw[as.character(i), as.character(j)] == 1) {
                scite.tree[i,j] = 1
            }
        }
    }
    colnames(scite.tree) = names[[name]]
    rownames(scite.tree) = names[[name]]
    print(scite.tree)

    fontsize = rep(44, ncol(scite.tree))
    names(fontsize) = names[[name]]

    graph = graph.adjacency(scite.tree)
    nel = igraph.to.graphNEL(graph)
    dev.new()
    plot(nel, nodeAttrs = list(fontsize = fontsize))
    title(paste0('\n', name))
    dev.copy2pdf(file = paste0('plot/', name, '.pdf'))
    dev.off()

}
