load('RData/dataset.single.cells.medium.RData')

dataset = dataset.single.cells.medium['50', , drop=F]
mini_dataset = dataset

epos = c(0.003, 0.004, 0.005, 0.006, 0.007)
eneg = c(0.03, 0.04, 0.05, 0.06, 0.07) 

for (i in 1:ncol(dataset)) {
	this_dataset = dataset[[1,i]][[2]]
	new_list_of_dataset = list()
	for (j in 1:length(epos)) {
		this_dataset$epos = epos[[j]]
		this_dataset$eneg = eneg[[j]]
		new_list_of_dataset[[as.character(j)]] = this_dataset
	}
	mini_dataset[[1,i]] = new_list_of_dataset
}

save('mini_dataset', file='RData/mini_dataset.R')
source('../generate.scite.input.Rdata')
create.scite.input(mini_dataset, 'single', 'mini', scite.sd)


### after recon

source('../reconstruct.scite.import.R')
source('../reconstruct.run.R')

load('RData/result.single.cells.mini.RData')

library(igraph)
library(sna)
library(Rgraphviz)

#### merge tronco results with scite
experiments.single.cells.mini.scite = import.scite.output(result.single.cells.mini, 'single', 'mini')
save(experiments.single.cells.mini.scite, file = 'RData/experiments.single.cells.mini.scite.RData')

### statistics

source('../statistics.plot.R')
source('../statistics.compute.R')
source('../giulio.plot.R')

experiments.single.cells.mini.scite.stats = get.stats(experiments.single.cells.mini.scite)
save(experiments.single.cells.mini.scite.stats, file="RData/experiments.single.cells.mini.scite.stats.RData")
giulio.plot(experiments.single.cells.mini.scite.stats, 'single', 'mini')

library(ggplot2)
library(Rmisc)

source('../giulio.plot.R')
samples = c(50)
e = new.env()

for (type in c('accuracy', 'hamming_distance', 'sensitivity', 'specificity')) {
    for (branching in c('mini')) {
        load(paste0('RData/results.values.', type, '.', branching, '.RData'), envir = e)
        results.values = e$results.values
        load(paste0('RData/results.', type, '.', branching, '.RData'), envir = e)
        results = e$results


        # ALL
        plotlist = list()
        plotlist.jitter = list()
        plot.id = 1

        for (sample in samples) {
            p = dotplotter(results.values, sample, c('capri_bic',
                'caprese_no.reg',
                'edmonds_pmi.no.reg',
                'gabow_pmi.no.reg',
                'chowliu_loglik', 'prim_no.reg',
                'scite_no.reg'), 
                branching,
                type,
                paste('SAMPLE SIZE = ', sample))

            plotlist[[plot.id]] = p
            
            m = jitterplotter(results.values, sample, c('capri_bic',
                'caprese_no.reg',
                'edmonds_pmi.no.reg',
                'gabow_pmi.no.reg',
                'chowliu_loglik', 'prim_no.reg',
                'scite_no.reg'), 
                branching,
                type,
                paste('SAMPLE SIZE = ', sample))
            plotlist.jitter[[plot.id]] = m
            plot.id = plot.id + 1
        }
        pdf(paste('plot_reduced/mt', type, branching, '.pdf', sep='_'), height = 5, width = 15)
        multiplot(plotlist = plotlist)
        dev.off()

        pdf(paste('plot_reduced/mt_jitter', type, branching, '.pdf', sep='_'), height = 5, width = 15)
        multiplot(plotlist = plotlist.jitter)
        dev.off()



    }
}
