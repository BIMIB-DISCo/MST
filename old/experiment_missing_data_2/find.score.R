library(bnlearn)
library(igraph)
library(Rgraphviz)
library(sna)

load('RData/experiment.missing.data.RData')

find.score = NULL
find.score$capri$id = 0
find.score$capri$score = NA
find.score$caprese$id = 0
find.score$caprese$score = NA
find.score$prim$id = 0
find.score$prim$score = NA
find.score$chowliu$id = 0
find.score$chowliu$score = NA
find.score$edmonds$id = 0
find.score$edmonds$score = NA
find.score$gabow$id = 0
find.score$gabow$score = NA
find.score$scite$score = NA

dataset = experiment.missing.data[[1,1]][[1]]
node.names = colnames(dataset$dataset)

plot_graph = function(adj, title, file) {
    colnames(adj) = node.names
    rownames(adj) = node.names
    fontsize = rep(20, ncol(adj))
    names(fontsize) = colnames(adj)
    graph = graph.adjacency(adj)
    nel = igraph.to.graphNEL(graph)
    dev.new()
    plot(nel, nodeAttrs = list(fontsize = fontsize))
    title(title)
    dev.copy2pdf(file = paste0('plot/', file, '.pdf'))
    dev.off()
}



for(sample in 1:nrow(experiment.missing.data)) {
    for (exp in 1:ncol(experiment.missing.data)) {
        exec = experiment.missing.data[[sample, exp]]
        for (noise in 1:length(exec)) {
            experiment = exec[[noise]]
            # capri
            if (is.na(find.score$capri$score) || experiment$reconstruction$capri$score > find.score$capri$score) {
                find.score$capri$score = experiment$reconstruction$capri$score
                find.score$capri$id = exp
                title = paste('\nCAPRI score:', round(find.score$capri$score, 3), ' id: ', find.score$capri$id)
                plot_graph(experiment$reconstruction$capri$bic.adj, title, 'capri')
            } else if (experiment$reconstruction$capri$score == find.score$capri$score) {
                find.score$capri$id = c(find.score$capri$id, exp)
            }

            # caprese
            if (is.na(find.score$caprese$score) || experiment$reconstruction$caprese$score > find.score$caprese$score) {
                find.score$caprese$score = experiment$reconstruction$caprese$score
                find.score$caprese$id = exp
                title = paste('\nCAPRESE score:', round(find.score$caprese$score, 3), ' id: ', find.score$caprese$id)
                plot_graph(experiment$reconstruction$capre$adj.caprese, title, 'caprese')
            } else if (experiment$reconstruction$caprese$score == find.score$caprese$score) {
                find.score$caprese$id = c(find.score$caprese$id, exp)
            }

            # prim
            if (is.na(find.score$prim$score) || experiment$reconstruction$prim$score > find.score$prim$score) {
                find.score$prim$score = experiment$reconstruction$prim$score
                find.score$prim$id = exp
                title = paste('\nPRIM score:', round(find.score$prim$score, 3), ' id: ', find.score$prim$id)
                plot_graph(experiment$reconstruction$prim$no.reg.adj, title, 'prim')
            } else if (experiment$reconstruction$prim$score == find.score$prim$score) {
                find.score$prim$id = c(find.score$prim$id, exp)
            }

            # chowliu
            if (is.na(find.score$chowliu$score) || experiment$reconstruction$chowliu$score > find.score$chowliu$score) {
                find.score$chowliu$score = experiment$reconstruction$chowliu$score
                find.score$chowliu$id = exp
                title = paste('\nCHOW LIU score:', round(find.score$chowliu$score, 3), ' id: ', find.score$chowliu$id)
                plot_graph(experiment$reconstruction$chowliu$loglik.adj, title, 'chowliu')
            } else if (experiment$reconstruction$chowliu$score == find.score$chowliu$score) {
                find.score$chowliu$id = c(find.score$chowliu$id, exp)
            }

            # edmonds
            if (is.na(find.score$edmonds$score) || experiment$reconstruction$edmonds$score > find.score$edmonds$score) {
                find.score$edmonds$score = experiment$reconstruction$edmonds$score
                find.score$edmonds$id = exp
                title = paste('\nEDMONDS score:', round(find.score$edmonds$score, 3), ' id: ', find.score$edmonds$id)
                plot_graph(experiment$reconstruction$edmonds$pmi.no.reg.adj, title, 'edmonds')
            } else if (experiment$reconstruction$edmonds$score == find.score$edmonds$score) {
                find.score$edmonds$id = c(find.score$edmonds$id, exp)
            }

            # gabow
            if (is.na(find.score$gabow$score) || experiment$reconstruction$gabow$score > find.score$gabow$score) {
                find.score$gabow$score = experiment$reconstruction$gabow$score
                find.score$gabow$id = exp
                title = paste('\nGABOW score:', round(find.score$gabow$score, 3), ' id: ', find.score$gabow$id)
                plot_graph(experiment$reconstruction$gabow$pmi.no.reg.adj, title, 'gabow')
            } else if (experiment$reconstruction$gabow$score == find.score$gabow$score) {
                find.score$gabow$id = c(find.score$gabow$id, exp)
            }
        }
    }
}

filename = paste0('scite_output/datasets/missing/single/58_1_1_ml0')
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

dataset = experiment.missing.data[[1,1]][[1]]
dataset = dataset$reconstructions$capri$bic.adj
colnames(scite.tree) = colnames(dataset)
rownames(scite.tree) = colnames(dataset)

load('RData/scite.RData')
colnames(scite) = colnames(dataset)

# create the igraph structure
net = empty.graph(colnames(scite), num = 4)
categoric.dataset = data.frame(apply(scite, 2, factor))
for (name in colnames(categoric.dataset)) {
    levels(categoric.dataset[[name]]) = c(0,1,3)
}
amat(net[[1]]) = scite.tree
find.score$scite$score = bnlearn::score(net[[1]], categoric.dataset, type='loglik') 
title = paste('\nSCITE score:', round(find.score$scite$score, 3))
plot_graph(scite.tree, title, 'scite')