##############################################################################
###
### MST
###
### Plot Random Main
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

samples = c(100, 500, 1000, 5000, 10000)
source('../giulio.plot.R')
e = new.env()

for (type in c('sensitivity', 'specificity', 'elasped_time')) {
    for (branching in c('random_20')) {
        load(paste0('RData/results.values.', type, '.', branching, '.RData'),
             envir = e)
        results.values = e$results.values
        load(paste0('RData/results.', type, '.', branching, '.RData'),
             envir = e)
        results = e$results

        ## EDMONDS
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        for (sample in samples) {
            p = dotplotter(results.values,
                           sample,
                           c('edmonds_entropy.no.reg',
                             'edmonds_pmi.no.reg',
                             'edmonds_cpmi.no.reg'),
                           'large_20',
                           type,
                           paste('SAMPLE SIZE = ', sample),
                           noise = c(1),
                           external.epos = c(0.005),
                           external.eneg = c(0.050))

            plotlist[[plot.id]] = p

            ##m = medianplotter(results,
            ##    sample,
            ##    c('edmonds_entropy.no.reg',
            ##    'edmonds_pmi.no.reg',
            ##    'edmonds_cpmi.no.reg'),
            ##    branching,
            ##    paste('SAMPLE SIZE = ', sample))
            ##plotlist.median[[plot.id]] = m
            
            plot.id = plot.id + 1
        }
        pdf(paste('plot_reduced/edmonds', type, branching, '.pdf',
                  sep = '_'),
            height = 10,
            width = 9)
        multiplot(plotlist = plotlist)
        dev.off()

        ##pdf(paste('plot_reduced/edmonds', branching, '_median.pdf', sep = '_'), height = 14, width = 11)
        ##multiplot(plotlist = plotlist.median)
        ##dev.off()

        ## CHOW LIU
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        for (sample in samples) {
            p = dotplotter(results.values,
                           sample,
                           c('chowliu_loglik'),
                           'large_20',
                           type,
                           paste('SAMPLE SIZE = ', sample),
                           noise = c(1),
                           external.epos = c(0.005),
                           external.eneg = c(0.050))

            plotlist[[plot.id]] = p

            ##m = medianplotter(results,
            ##    sample,
            ##    c('chowliu_loglik'),
            ##    branching,
            ##    paste('SAMPLE SIZE = ', sample))
            ##plotlist.median[[plot.id]] = m
            plot.id = plot.id + 1
        }
        pdf(paste('plot_reduced/chowliu', type, branching, '.pdf',
                  sep = '_'),
            height = 10,
            width = 9)
        multiplot(plotlist = plotlist)
        dev.off()

        ##pdf(paste('plot_reduced/prim_chowliu', branching, '_median.pdf', sep = '_'), height = 14, width = 11)
        ##multiplot(plotlist = plotlist.median)
        ##dev.off()


        ## ALL
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        for (sample in samples) {
            p = dotplotter(results.values,
                           sample,
                           c('edmonds_pmi.no.reg',
                             'chowliu_loglik'),
                           'large_20',
                           type,
                           paste('SAMPLE SIZE = ', sample),
                           noise = c(1),
                           external.epos = c(0.005),
                           external.eneg = c(0.050))

            plotlist[[plot.id]] = p

            ##m = medianplotter(results,
            ##    sample,
            ##    c('edmonds_entropy.no.reg',
            ##    'edmonds_pmi.no.reg',
            ##    'edmonds_cpmi.no.reg'),
            ##    branching,
            ##    paste('SAMPLE SIZE = ', sample))
            ##plotlist.median[[plot.id]] = m
            
            plot.id = plot.id + 1
        }
        pdf(paste('plot_reduced/all', type, branching, '.pdf', sep = '_'), height = 10, width = 9)
        multiplot(plotlist = plotlist)
        dev.off()
    }
}

### end of file -- plot.random.main.reduced.R
