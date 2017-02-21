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
    model = change.color(model, type='variant', new.color='lightblue3')
    oncoprint(model, file = paste0('plot/', file, '_data.pdf'), title = title)
    tronco.plot(model, 
        edge.cex = 1.5,          # scale edge size
        legend.cex = .5,         # scale legend size
        scale.nodes = .6,        # scale node size
        confidence = c('tp', 'pr', 'hg'), # display p-values for these statistics 
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
                title = paste('\nCAPRI score:', round(find.score$capri$score, 3), ' id: ', find.score$capri$id)
                plot_graph(experiment$reconstruction$capri$model, title, 'capri')
                capri_best_model = experiment$reconstruction$capri$model
                save(capri_best_model,file="capri_best_model.RData")
            } else if (experiment$reconstruction$capri$score == find.score$capri$score) {
                find.score$capri$id = c(find.score$capri$id, exp)
            }

            # caprese
            if (is.na(find.score$caprese$score) || experiment$reconstruction$caprese$score > find.score$caprese$score) {
                find.score$caprese$score = experiment$reconstruction$caprese$score
                find.score$caprese$id = exp
                title = paste('\nCAPRESE score:', round(find.score$caprese$score, 3), ' id: ', find.score$caprese$id)
                plot_graph(experiment$reconstruction$capre$model, title, 'caprese')
                caprese_best_model = experiment$reconstruction$caprese$model
                save(caprese_best_model,file="caprese_best_model.RData")
            } else if (experiment$reconstruction$caprese$score == find.score$caprese$score) {
                find.score$caprese$id = c(find.score$caprese$id, exp)
            }

            # prim
            if (is.na(find.score$prim$score) || experiment$reconstruction$prim$score > find.score$prim$score) {
                find.score$prim$score = experiment$reconstruction$prim$score
                find.score$prim$id = exp
                title = paste('\nPRIM score:', round(find.score$prim$score, 3), ' id: ', find.score$prim$id)
                plot_graph(experiment$reconstruction$prim$model, title, 'prim')
                prim_best_model = experiment$reconstruction$prim$model
                save(prim_best_model,file="prim_best_model.RData")
            } else if (experiment$reconstruction$prim$score == find.score$prim$score) {
                find.score$prim$id = c(find.score$prim$id, exp)
            }

            # chowliu
            if (is.na(find.score$chowliu$score) || experiment$reconstruction$chowliu$score > find.score$chowliu$score) {
                find.score$chowliu$score = experiment$reconstruction$chowliu$score
                find.score$chowliu$id = exp
                title = paste('\nCHOW LIU score:', round(find.score$chowliu$score, 3), ' id: ', find.score$chowliu$id)
                plot_graph(experiment$reconstruction$chowliu$model, title, 'chowliu')
                chowliu_best_model = experiment$reconstruction$chowliu$model
                save(chowliu_best_model,file="chowliu_best_model.RData")
            } else if (experiment$reconstruction$chowliu$score == find.score$chowliu$score) {
                find.score$chowliu$id = c(find.score$chowliu$id, exp)
            }

            # edmonds
            if (is.na(find.score$edmonds$score) || experiment$reconstruction$edmonds$score > find.score$edmonds$score) {
                find.score$edmonds$score = experiment$reconstruction$edmonds$score
                find.score$edmonds$id = exp
                title = paste('\nEDMONDS score:', round(find.score$edmonds$score, 3), ' id: ', find.score$edmonds$id)
                plot_graph(experiment$reconstruction$edmonds$model, title, 'edmonds')
                edmonds_best_model = experiment$reconstruction$edmonds$model
                save(edmonds_best_model,file="edmonds_best_model.RData")
            } else if (experiment$reconstruction$edmonds$score == find.score$edmonds$score) {
                find.score$edmonds$id = c(find.score$edmonds$id, exp)
            }

            # gabow
            if (is.na(find.score$gabow$score) || experiment$reconstruction$gabow$score > find.score$gabow$score) {
                find.score$gabow$score = experiment$reconstruction$gabow$score
                find.score$gabow$id = exp
                title = paste('\nGABOW score:', round(find.score$gabow$score, 3), ' id: ', find.score$gabow$id)
                plot_graph(experiment$reconstruction$gabow$model, title, 'gabow')
                gabow_best_model = experiment$reconstruction$gabow$model
                save(gabow_best_model,file="gabow_best_model.RData")
            } else if (experiment$reconstruction$gabow$score == find.score$gabow$score) {
                find.score$gabow$id = c(find.score$gabow$id, exp)
            }
        }
    }
}
