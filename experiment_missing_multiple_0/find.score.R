library(bnlearn)
library(igraph)
library(Rgraphviz)
library(sna)
library(TRONCO)

load('RData/experiment.missing.data.RData')
source('../reconstruct.run.R')


#plot_graph = function(model, title, file) {
#    model = change.color(model, type='variant', new.color='lightgreen')
#    oncoprint(model, file = paste0('plot/', file, '_data.pdf'), title = title)
#    tronco.plot(model, 
#        edge.cex = 1.5,          # scale edge size
#        legend.cex = .5,         # scale legend size
#        scale.nodes = .6,        # scale node size
#        confidence = c('tp', 'pr', 'hg', 'npb'), # display p-values for these statistics 
#        disconnected = FALSE,        # do not display nodes without incoming/outgoing edges
#        height.logic = .3,       # scale logical connectives
#        file = paste0('plot/', file, '.pdf'),
#        title = title)
#}

results.experiment.missing.data = experiment.missing.data

for(sample in 1:nrow(experiment.missing.data)) {
    for (exp in 1:ncol(experiment.missing.data)) {
        exec = experiment.missing.data[[sample, exp]]
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
        this.find.score = find.score
        for (imputated in 1:length(exec)) {
            experiment = exec[[imputated]]

            # create a TRONCO object of the dataset
            data = import.genotypes(experiment$dataset)

            # create the igraph structure
            net = empty.graph(colnames(as.genotypes(data)), num = 6)
            categoric.dataset = data.frame(apply(as.genotypes(data), 2, factor))
            for (name in colnames(categoric.dataset)) {
                levels(categoric.dataset[[name]]) = c(0,1)
            }

            # capri
            amat(net[[1]]) = experiment$reconstructions$capri$adj
            experiment$reconstruction$capri$score = bnlearn::score(net[[1]], categoric.dataset, type='loglik')
            if (is.na(this.find.score$capri$score) || experiment$reconstruction$capri$score > this.find.score$capri$score) {
                this.find.score$capri$score = experiment$reconstruction$capri$score
                this.find.score$capri$id = imputated
            }

            # caprese
            amat(net[[2]]) = experiment$reconstructions$caprese$adj
            experiment$reconstruction$caprese$score = bnlearn::score(net[[2]], categoric.dataset, type='loglik')
            if (is.na(this.find.score$caprese$score) || experiment$reconstruction$caprese$score > this.find.score$caprese$score) {
                this.find.score$caprese$score = experiment$reconstruction$caprese$score
                this.find.score$caprese$id = imputated
            }

            # prim
            amat(net[[3]]) = experiment$reconstructions$prim$adj
            experiment$reconstruction$prim$score = bnlearn::score(net[[3]], categoric.dataset, type='loglik')
            if (is.na(this.find.score$prim$score) || experiment$reconstruction$prim$score > this.find.score$prim$score) {
                this.find.score$prim$score = experiment$reconstruction$prim$score
                this.find.score$prim$id = imputated
            }

            # chowliu
            amat(net[[4]]) = experiment$reconstructions$chowliu$loglik.adj
            experiment$reconstruction$chowliu$score = bnlearn::score(net[[4]], categoric.dataset, type='loglik')
            if (is.na(this.find.score$chowliu$score) || experiment$reconstruction$chowliu$score > this.find.score$chowliu$score) {
                this.find.score$chowliu$score = experiment$reconstruction$chowliu$score
                this.find.score$chowliu$id = imputated
            }

            # edmonds
            amat(net[[5]]) = experiment$reconstructions$edmonds$adj
            experiment$reconstruction$edmonds$score = bnlearn::score(net[[5]], categoric.dataset, type='loglik')
            if (is.na(this.find.score$edmonds$score) || experiment$reconstruction$edmonds$score > this.find.score$edmonds$score) {
                this.find.score$edmonds$score = experiment$reconstruction$edmonds$score
                this.find.score$edmonds$id = imputated
            }

            # gabow
            amat(net[[6]]) = experiment$reconstructions$gabow$pmi.no.reg.adj
            experiment$reconstruction$gabow$score = bnlearn::score(net[[6]], categoric.dataset, type='loglik')
            if (is.na(this.find.score$gabow$score) || experiment$reconstruction$gabow$score > this.find.score$gabow$score) {
                this.find.score$gabow$score = experiment$reconstruction$gabow$score
                this.find.score$gabow$id = imputated
            }
        }
        results.experiment.missing.data[[sample, exp]] = this.find.score
    }
}

best.score = matrix(list(), ncol = ncol(results.experiment.missing.data), nrow = nrow(results.experiment.missing.data))
colnames(best.score) =  colnames(results.experiment.missing.data)
rownames(best.score) =  rownames(results.experiment.missing.data)

for(sample in 1:nrow(experiment.missing.data)) {
    for (exp in 1:ncol(experiment.missing.data)) {
        this.result = NULL
        for (algo in names(results.experiment.missing.data[[sample, exp]])) {
            print(algo)
            id = results.experiment.missing.data[[sample, exp]][[algo]][['id']]
            print(id)
            this.result[[algo]][['dataset']] = experiment.missing.data[[sample, exp]][[id]][['dataset']]
            this.result[[algo]][['true_tree']] = experiment.missing.data[[sample, exp]][[id]][['true_tree']]
            this.result[[algo]][['epos']] = experiment.missing.data[[sample, exp]][[id]][['epos']]
            this.result[[algo]][['eneg']] = experiment.missing.data[[sample, exp]][[id]][['eneg']]
            this.result[[algo]][['reconstructions']] = experiment.missing.data[[sample, exp]][[id]][['reconstructions']][[algo]]
        }
        best.score[[sample, exp]] = this.result
    }
}


prefix = paste0('scite_output/datasets/missing/multiple/')
for(sample in 1:nrow(best.score)) {
    for (exp in 1:ncol(best.score)) {
        this.best.score = best.score[[sample, exp]]
        name = paste0('20_', exp, '_', sample, '_ml0')
        complete.filename = paste0(prefix, name)
        read = readLines(paste0(complete.filename, '.gv'))
        read = gsub(' ', '', read)
        read = gsub(';', '', read)
        write(paste0(read, collapse = '\n'), file = paste0(complete.filename, '.correct.gv'))
        readdot = read.dot(paste0(complete.filename, '.correct.gv'))
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
        true.tree = this.best.score[[1]][['true_tree']]
        epos = this.best.score[[1]][['epos']]
        eneg = this.best.score[[1]][['eneg']]
        scite.result = NULL
        scite.result$true_tree = true.tree
        scite.result$epos = epos
        scite.result$eneg = eneg
        scite.result$reconstructions$result = getStats(true.tree, scite.tree)
        scite.result$reconstructions$adj = scite.tree
        this.best.score$scite = scite.result
        best.score[[sample, exp]] = this.best.score
    }
}

save(best.score, file='RData/best.score.RData')



create.data.frame <- function(dataset) {
    results = data.frame(x = NULL, stringsAsFactors = FALSE)

    for(sample in 1:nrow(best.score)) {
        for (exp in 1:ncol(best.score)) {
            this.exp = best.score[[sample, exp]]
            for (algo in names(this.exp)) {


                this.exp.algo = this.exp[[algo]]
                results = rbind(results, c(sample,
                            exp,
                            this.exp.algo[['reconstructions']][['result']][['sensitivity']],
                            this.exp.algo[['reconstructions']][['result']][['specificity']],
                            algo,
                            paste0('MD_',sample)), stringsAsFactors = FALSE)

            }
        }
    }


    colnames(results) = c('sample', 'exp', 'sensitivity', 'specificity', 'algorithm', 'source')

    results$sample = as.factor(results$sample)
    results$exp = as.factor(results$exp)
    results$sensitivity = as.numeric(results$sensitivity)
    results$specificity = as.numeric(results$specificity)
    results$algorithm = as.factor(results$algorithm)
    results$source = as.factor(results$source)

    save(results, file=paste0('RData/results.missing.data.RData'))
    return(results)
}

results.missing.data = create.data.frame(best.score)

load('RData/experiments.multiple.biopses.medium.scite.stats.RData')

reg = NULL
reg$capri = 'bic'
reg$gabow = 'pmi.no.reg'
reg$caprese = 'no.reg'
reg$prim = 'no.reg'
reg$edmonds = 'pmi.no.reg'
reg$chowliu = 'loglik'
reg$scite = 'no.reg'

results.original.data = data.frame(x = NULL, stringsAsFactors = FALSE)

for (algo in c('caprese', 'capri', 'prim', 'edmonds', 'gabow', 'prim', 'chowliu', 'scite')) {
    this.reg = reg[[algo]]
    sensitivity = experiments.single.cells.medium.scite.stats[['sensitivity']][[algo]][[this.reg]][[1,4]][['values']][1:10]
    specificity = experiments.single.cells.medium.scite.stats[['specificity']][[algo]][[this.reg]][[1,4]][['values']][1:10]
    for (i in 1:10) {
        results.original.data = rbind(results.original.data, c(0,
            i,
            sensitivity[[i]],
            specificity[[i]],
            algo,
            'EXP2'), stringsAsFactors = FALSE)
    }
    
}

colnames(results.original.data) = c('sample', 'exp', 'sensitivity', 'specificity', 'algorithm', 'source')

results.original.data$sample = as.factor(results.original.data$sample)
results.original.data$exp = as.factor(results.original.data$exp)
results.original.data$sensitivity = as.numeric(results.original.data$sensitivity)
results.original.data$specificity = as.numeric(results.original.data$specificity)
results.original.data$algorithm = as.factor(results.original.data$algorithm)
results.original.data$source = as.factor(results.original.data$source)

save(results.original.data, file='RData/results.original.data.RData')


complete.results = rbind(results.original.data, results.missing.data)
save(complete.results, file = 'RData/complete.results.RData')

experiment.names = c('EXP2' = 'Original dataset',
    'MD_1' = 'Missing data 10%',                           
    'MD_2' = 'Missing data 20%',
    'MD_3' = 'Missing data 30%',
    'MD_4' = 'Missing data 40%')

experiment.palette = c('EXP2' = 'red',
    'MD_1' = '#fa9fb5',                           
    'MD_2' = '#f768a1',
    'MD_3' = '#dd3497',
    'MD_4' = '#ae017e')

description = c('sensitivity' = 'Sensitivity',
    'specificity' = 'Specificity')

library(ggplot2)

p = ggplot(complete.results, aes(x = algorithm, y = sensitivity, fill = source)) +
    scale_fill_manual(values = experiment.palette, labels = experiment.names) +
    ylab(paste0('Sensitivity (n = 11)')) +
    geom_boxplot(outlier.size = 0) +
    xlab('epos = 0.0; eneg = 0.0')

dev.new(width = 16, height = 5)
print(p)
dev.copy2pdf(file = 'plot/sensitivity.pdf')
dev.off()


p = ggplot(complete.results, aes(x = algorithm, y = specificity, fill = source)) +
    scale_fill_manual(values = experiment.palette, labels = experiment.names) +
    ylab(paste0('Specificity (n = 11)')) +
    geom_boxplot(outlier.size = 0) +
    xlab('epos = 0.0; eneg = 0.0')

dev.new(width = 16, height = 5)
print(p)
dev.copy2pdf(file = 'plot/specificity.pdf')
dev.off()


