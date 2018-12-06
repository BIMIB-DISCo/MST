##############################################################################
###
### MST
###
### Experiment 4 Multiple
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

load('RData/dataset.multiple.biopses.medium.RData')

dataset = dataset.multiple.biopses.medium['20', , drop = F]
mini_dataset = dataset

epos = c(0.03, 0.04, 0.05, 0.06, 0.07)
eneg = c(0.03, 0.04, 0.05, 0.06, 0.07) 

for (i in 1:ncol(dataset)) {
    this_dataset = dataset[[1, i]][[2]]
    new_list_of_dataset = list()
    for (j in 1:length(epos)) {
        this_dataset$epos = epos[[j]]
        this_dataset$eneg = eneg[[j]]
        new_list_of_dataset[[as.character(j)]] = this_dataset
    }
    mini_dataset[[1,i]] = new_list_of_dataset
}

save('mini_dataset', file = 'RData/mini_dataset.RData')
source('../generate.scite.input.R')
create.scite.input(mini_dataset, 'multiple', 'mini', scite.sd)


### After recon

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/result.multiple.biopses.mini.RData')

library(igraph)
library(sna)
library(Rgraphviz)


### Merge tronco results with scite

experiments.multiple.biopses.mini.scite =
    import.scite.output(result.multiple.biopses.mini, 'multiple', 'mini')
save(experiments.multiple.biopses.mini.scite,
     file = 'RData/experiments.multiple.biopses.mini.scite.RData')


### Statistics

source('../statistics.plot.R')
source('../statistics.compute.R')
source('../giulio.plot.R')

experiments.multiple.biopses.mini.scite.stats =
    get.stats(experiments.multiple.biopses.mini.scite)
save(experiments.multiple.biopses.mini.scite.stats,
     file = "RData/experiments.multiple.biopses.mini.scite.stats.RData")
giulio.plot(experiments.multiple.biopses.mini.scite.stats, 'multiple', 'mini')

library(ggplot2)
library(Rmisc)

source('../giulio.plot.R')
samples = c(20)
e = new.env()

for (type in c('accuracy', 'hamming_distance', 'sensitivity', 'specificity')) {
    for (branching in c('mini')) {
        load(paste0('RData/results.values.', type, '.', branching, '.RData'), envir = e)
        results.values = e$results.values
        load(paste0('RData/results.', type, '.', branching, '.RData'), envir = e)
        results = e$results


        ## ALL
        plotlist = list()
        plotlist.jitter = list()
        plot.id = 1

        for (sample in samples) {
            p = dotplotter(results.values,
                           sample,
                           c('capri_bic', 
                             'caprese_no.reg', 
                             'edmonds_pmi.no.reg', 
                             'gabow_pmi.no.reg', 
                             'chowliu_loglik',
                             'prim_no.reg', 
                             'scite_no.reg'), 
                           branching, 
                           type, 
                           paste('SAMPLE SIZE = ', sample), 
                           sample.type = 'multiple')

            plotlist[[plot.id]] = p
            
            m = jitterplotter(results.values,
                              sample,
                              c('capri_bic', 
                                'caprese_no.reg', 
                                'edmonds_pmi.no.reg', 
                                'gabow_pmi.no.reg', 
                                'chowliu_loglik',
                                'prim_no.reg', 
                                'scite_no.reg'), 
                              branching, 
                              type, 
                              paste('SAMPLE SIZE = ', sample), 
                              sample.type = 'multiple')
            plotlist.jitter[[plot.id]] = m
            plot.id = plot.id + 1
        }
        pdf(paste('plot_reduced/mt', type, branching, '.pdf',
                  sep = '_'),
            height = 5,
            width = 15)
        multiplot(plotlist = plotlist)
        dev.off()

        pdf(paste('plot_reduced/mt_jitter', type, branching, '.pdf',
                  sep = '_'),
            height = 5,
            width = 15)
        multiplot(plotlist = plotlist.jitter)
        dev.off()
    }
}

### end of file -- experiment_4.R
