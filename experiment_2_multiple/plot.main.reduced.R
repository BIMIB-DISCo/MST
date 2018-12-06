##############################################################################
###
### MST
###
### Plot Main Reduced
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

library(ggplot2)
library(Rmisc)

samples = c(5, 10, 20)
source('../giulio.plot.R')
e = new.env()

for (type in c('accuracy', 'hamming_distance', 'sensitivity', 'specificity')) {
    for (branching in c('low', 'medium', 'high')) {
        load(paste0('RData/results.values.', type, '.', branching, '.RData'),
             envir = e)
        results.values = e$results.values
        load(paste0('RData/results.', type, '.', branching, '.RData'),
             envir = e)
        results = e$results

        ## ALL
        plotlist = list()
        plotlist.median = list()
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
                           noise = c(1,2,5))

            plotlist[[plot.id]] = p
            
            m = medianplotter(results,
                              sample,
                              c('edmonds_entropy.no.reg',
                                'edmonds_pmi.no.reg',
                                'edmonds_cpmi.no.reg'), 
                              branching,
                              paste('SAMPLE SIZE = ', sample))
            plotlist.median[[plot.id]] = m
            plot.id = plot.id + 1
        }
        pdf(paste('plot_reduced/all', type, branching, '.pdf',
                  sep = '_'),
            height = 10,
            width = 9)
        multiplot(plotlist = plotlist)
        dev.off()
    }
}



samples = c(10)
e = new.env()

for (type in c('accuracy', 'hamming_distance', 'sensitivity', 'specificity')) {
    for (branching in c('medium')) {
        load(paste0('RData/results.values.', type, '.', branching, '.RData'),
             envir = e)
        results.values = e$results.values
        load(paste0('RData/results.', type, '.', branching, '.RData'),
             envir = e)
        results = e$results


        ## ALL
        plotlist = list()
        plotlist.jitter = list()
        plot.id = 1

        for (sample in samples) {
            p = dotplotter(results.values,
                           sample,
                           c('capri_bic',
                             'edmonds_pmi.no.reg',
                             'scite_no.reg'), 
                           branching,
                           type,
                           paste('SAMPLE SIZE = ', sample),
                           noise = c(2),
                           sample.type = 'multiple')

            plotlist[[plot.id]] = p
            
            m = jitterplotter(results.values,
                              sample,
                              c('capri_bic',
                                'edmonds_pmi.no.reg',
                                'scite_no.reg'), 
                              branching,
                              type,
                              paste('SAMPLE SIZE = ', sample),
                              noise = c(2),
                              sample.type = 'multiple')
            plotlist.jitter[[plot.id]] = m
            plot.id = plot.id + 1
        }
        pdf(paste('plot_reduced/mt', type, branching, '.pdf',
                  sep = '_'),
            height = 5,
            width = 5)
        multiplot(plotlist = plotlist)
        dev.off()

        pdf(paste('plot_reduced/mt_jitter', type, branching, '.pdf',
                  sep = '_'),
            height = 5,
            width = 5)
        multiplot(plotlist = plotlist.jitter)
        dev.off()
    }
}

### end of file -- plot.main.reduced.R
