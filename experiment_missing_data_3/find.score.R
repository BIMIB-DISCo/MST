library(bnlearn)
library(igraph)
library(Rgraphviz)
library(sna)
library(TRONCO)

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

dataset = experiment.missing.data[[1,1]][[1]]
node.names = colnames(dataset$dataset)

plot_graph = function(model, title, file) {
    model = change.color(model, type='variant', new.color='lightgreen')
    oncoprint(model, file = paste0('plot/', file, '_data.pdf'), title = title)
    tronco.plot(model, 
        edge.cex = 1.5,          # scale edge size
        legend.cex = .5,         # scale legend size
        scale.nodes = .6,        # scale node size
        confidence = c('tp', 'pr', 'hg', 'npb'), # display p-values for these statistics 
        disconnected = FALSE,        # do not display nodes without incoming/outgoing edges
        height.logic = .3,       # scale logical connectives
        file = paste0('plot/', file, '.pdf'),
        title = title)
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
            } else if (experiment$reconstruction$capri$score == find.score$capri$score) {
                find.score$capri$id = c(find.score$capri$id, exp)
            }

            # caprese
            if (is.na(find.score$caprese$score) || experiment$reconstruction$caprese$score > find.score$caprese$score) {
                find.score$caprese$score = experiment$reconstruction$caprese$score
                find.score$caprese$id = exp
            } else if (experiment$reconstruction$caprese$score == find.score$caprese$score) {
                find.score$caprese$id = c(find.score$caprese$id, exp)
            }

            # prim
            if (is.na(find.score$prim$score) || experiment$reconstruction$prim$score > find.score$prim$score) {
                find.score$prim$score = experiment$reconstruction$prim$score
                find.score$prim$id = exp
            } else if (experiment$reconstruction$prim$score == find.score$prim$score) {
                find.score$prim$id = c(find.score$prim$id, exp)
            }

            # chowliu
            if (is.na(find.score$chowliu$score) || experiment$reconstruction$chowliu$score > find.score$chowliu$score) {
                find.score$chowliu$score = experiment$reconstruction$chowliu$score
                find.score$chowliu$id = exp
            } else if (experiment$reconstruction$chowliu$score == find.score$chowliu$score) {
                find.score$chowliu$id = c(find.score$chowliu$id, exp)
            }

            # edmonds
            if (is.na(find.score$edmonds$score) || experiment$reconstruction$edmonds$score > find.score$edmonds$score) {
                find.score$edmonds$score = experiment$reconstruction$edmonds$score
                find.score$edmonds$id = exp
            } else if (experiment$reconstruction$edmonds$score == find.score$edmonds$score) {
                find.score$edmonds$id = c(find.score$edmonds$id, exp)
            }

            # gabow
            if (is.na(find.score$gabow$score) || experiment$reconstruction$gabow$score > find.score$gabow$score) {
                find.score$gabow$score = experiment$reconstruction$gabow$score
                find.score$gabow$id = exp
            } else if (experiment$reconstruction$gabow$score == find.score$gabow$score) {
                find.score$gabow$id = c(find.score$gabow$id, exp)
            }
        }
    }
}


for (name in ls(find.score)) {
    ids = find.score[[name]][['id']]
    score = find.score[[name]][['score']]
    for (id in ids) {
        print(id)
        print(score)
        model = experiment.missing.data[[1,id]][[1]][['reconstructions']][[name]][['model']]
        model.bootstrap = tronco.bootstrap(model)
        find.score[[name]][['model']] = model.bootstrap
    }
}

save(find.score, file = 'RData/find.score.RData')

for (name in ls(find.score)) {
    ids = find.score[[name]][['id']]
    score = find.score[[name]][['score']]
    for (id in ids) {
        model.bootstrap = find.score[[name]][['model']]
        title = paste('\n',name, 'score:', round(score, 3), ' id: ', id)
        plot_graph(model.bootstrap, title, paste0(name, '_', id))
    }
}


filename = paste0('scite_output/datasets/missing/single/47_1_1_ml0')
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

load('RData/find.score.RData')
dataset = experiment.missing.data[[1,1]][[1]]
dataset = dataset$dataset
colnames(scite.tree) = colnames(dataset)
rownames(scite.tree) = colnames(dataset)
graph = graph.adjacency(scite.tree)
graph = igraph.to.graphNEL(graph)
edge_names = edgeNames(graph)
eAttrs = list()
eAttrs$label = rep('', length(edge_names))
names(eAttrs$label) = edge_names
eAttrs$fontsize = rep(50, length(edge_names))
names(eAttrs$fontsize) = edge_names

# count edges
for (from in rownames(scite.tree)) {
    for (to in colnames(scite.tree)) {
        if (scite.tree[[from,to]] == 1) {
            label = ''
            edge.name = paste0(from,'~',to)
            for (algo in ls(find.score)) {
                model.bootstrap = find.score[[algo]][['model']]
                adj = as.adj.matrix(model.bootstrap)[[1]]
                colnames(adj) = colnames(dataset)
                rownames(adj) = colnames(dataset)
                if (adj[[from,to]] == 1) {
                    label = paste0(label, ' ', algo)
                }
            }
            eAttrs$label[[edge.name]] = label
        }
    }
}



save(scite.tree, file = 'RData/scite.tree')
load('RData/scite.RData')
colnames(scite) = colnames(dataset)

# create the igraph structure
net = empty.graph(colnames(scite), num = 4)
categoric.dataset = data.frame(apply(scite, 2, factor))
for (name in colnames(categoric.dataset)) {
    levels(categoric.dataset[[name]]) = c(0,1,3)
}
amat(net[[1]]) = scite.tree
scite.score = bnlearn::score(net[[1]], categoric.dataset, type='loglik') 
title = paste('SCITE score:', round(scite.score, 3))
fontsize = rep(30, length(colnames(dataset)))
nAttrs = list()
nAttrs$fontsize = fontsize
names(nAttrs$fontsize) = colnames(dataset)
dev.new()
plot(graph, edgeAttrs = eAttrs, nodeAttrs = nAttrs, main = title)
dev.copy2pdf(file = 'plot/scite.pdf')
#dev.off()
