# load library
library(igraph)
library(bnlearn)
library(Rgraphviz)
library(sna)
library(oncoNEM)
library(TRONCO)


names = NULL

load('data.RData')
complete = data
geno = as.genotypes(complete)
geno = keysToNames(complete, geno)
names$complete = colnames(geno)

load('data.merged.RData')
ind = data.merged
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

scite.bootstrap = matrix(0, ncol = ncol(geno), nrow = ncol(geno))
colnames(scite.bootstrap) = colnames(geno)
rownames(scite.bootstrap) = colnames(geno)

for (i in 1:100) {
    filename = paste0('scite_output/datasets/bootstrap_', i, '_ml0')
    read = readLines(paste0(filename, '.gv'))
    read = gsub(' ', '', read)
    read = gsub(';', '', read)
    write(paste0(read, collapse = '\n'), file = paste0(filename, '.correct.gv'))
    readdot = read.dot(paste0(filename, '.correct.gv'))
    graph = graph.adjacency(readdot)
    scite.tree.raw = get.adjacency(graph, sparse = FALSE)
    
    for (i in 1:nrow(scite.bootstrap)) {
        for(j in 1:ncol(scite.bootstrap)) {
            if (scite.tree.raw[as.character(i), as.character(j)] == 1) {
                scite.bootstrap[i,j] = scite.bootstrap[[i,j]] + 1
            }
        }
    }
    print(scite.bootstrap)



}

fontsize = rep(44, ncol(scite.tree))
names(fontsize) = colnames(geno)

graph = graph.adjacency(scite.tree)
nel = igraph.to.graphNEL(graph)


edge_names = edgeNames(nel)
eAttrs = list()

eAttrs$lwd = rep(1, length(edge_names))
names(eAttrs$lwd) = edge_names

eAttrs$label = rep('', length(edge_names))
names(eAttrs$label) = edge_names

eAttrs$fontsize = rep(44, length(edge_names))
names(eAttrs$fontsize) = edge_names

for (from in rownames(scite.tree)) {
    for(to in colnames(scite.tree)) {
        if (scite.tree[[from, to]] == 1) {
            value = scite.bootstrap[[from, to]]
            edge_name = paste0(from, '~', to)
            eAttrs$label[edge_name] = paste0('          ', as.character(round(value / 100, 2)))
            eAttrs$lwd[edge_name] = 0.8 + (value * 0.10)
        }
    }
}

dev.new()
plot(nel, nodeAttrs = list(fontsize = fontsize), edgeAttrs = eAttrs)
title('\n bootstrap')
dev.copy2pdf(file = paste0('plot/bootstrap.pdf'))
dev.off()



