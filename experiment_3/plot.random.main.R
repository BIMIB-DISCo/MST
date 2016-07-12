library(ggplot2)
library(Rmisc)

source('../giulio.plot.R')
load('RData/results.values.random_5.RData')
results.values.random_5 = results.values
load('RData/results.values.random_10.RData')
results.values.random_10 = results.values
load('RData/results.values.random_15.RData')
results.values.random_15 = results.values
load('RData/results.values.random_20.RData')
results.values.random_20 = results.values

load('RData/results.random_5.RData')
results.random_5 = results
load('RData/results.random_10.RData')
results.random_10 = results
load('RData/results.random_15.RData')
results.random_15 = results
load('RData/results.random_20.RData')
results.random_20 = results


for (branching in c('random_5', 'random_10', 'random_15', 'random_20')) {
    if (branching == 'random_5') {
        results.values = results.values.random_5
        results = results.random_5
    } else if (branching == 'random_10') {
        results.values = results.values.random_10
        results = results.random_10
    } else if (branching == 'random_15') {
        results.values = results.values.random_15
        results = results.random_15
    } else if (branching == 'random_20') {
        results.values = results.values.random_20
        results = results.random_20
    }


    # CAPRI CAPRESE
    plotlist = list()
    plotlist.median = list()
    plot.id = 1

    cat('capri caprese', branching, '\n')

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, 
            sample,
            c('capri_bic', 'capri_aic', 'capri_loglik', 'caprese_no.reg', 'scite_no.reg'),
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist[[plot.id]] = p

        m = medianplotter(results,
            sample,
            c('capri_bic', 'capri_aic', 'capri_loglik', 'caprese_no.reg'),
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist.median[[plot.id]] = m
        plot.id = plot.id + 1
    }
    pdf(paste('plot/capri_caprese', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()

    #pdf(paste('plot/capri_caprese', branching, '_median.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist.median)
    #dev.off()

    # GABOW
    plotlist = list()
    plotlist.median = list()
    plot.id = 1

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, sample, c('gabow_entropy.no.reg', 'gabow_pmi.no.reg', 'gabow_cpmi.no.reg',
            'gabow_mi.no.reg', 'scite_no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist[[plot.id]] = p

        m = medianplotter(results,
            sample,c('gabow_entropy.no.reg', 'gabow_pmi.no.reg', 'gabow_cpmi.no.reg',
            'gabow_mi.no.reg', 'gabow_no_rising_no.raising.entropy.no.reg',
            'gabow_no_rising_no.raising.pmi.no.reg', 'gabow_no_rising_no.raising.cpmi.no.reg',
            'gabow_no_rising_no.raising.mi.no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist.median[[plot.id]] = m
        plot.id = plot.id + 1
    }
    pdf(paste('plot/gabow', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()

    #pdf(paste('plot/gabow', branching, '_median.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist.median)
    #dev.off()

    # EDMONDS
    plotlist = list()
    plotlist.median = list()
    plot.id = 1

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, 
            sample, 
            c('edmonds_entropy.no.reg',
            'edmonds_pmi.no.reg',
            'edmonds_cpmi.no.reg',
            'scite_no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))
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
    pdf(paste('plot/edmonds', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()

    #pdf(paste('plot/edmonds', branching, '_median.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist.median)
    #dev.off()

    # SCITE
    #plotlist = list()
    #plotlist.median = list()
    #plot.id = 1
    #
    #for (sample in c(10, 25, 50, 75, 100)) {
    #    p = dotplotter(results.values, sample, c('scite_no.reg'), 
    #        branching,
    #        paste('SAMPLE SIZE = ', sample))
    #    plotlist[[plot.id]] = p
    #    
    #    m = medianplotter(results,
    #        sample,
    #        c('scite_no.reg'), 
    #        branching,
    #        paste('SAMPLE SIZE = ', sample))
    #    plotlist.median[[plot.id]] = m
    #    plot.id = plot.id + 1
    #}
    #pdf(paste('plot/scite', branching, '.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist)
    #dev.off()

    #pdf(paste('plot/scite', branching, '_median.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist.median)
    #dev.off()

    # PRIM CHOW LIU
    plotlist = list()
    plotlist.median = list()
    plot.id = 1

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, sample, c('chowliu_loglik', 'prim_no.reg', 'scite_no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist[[plot.id]] = p

        m = medianplotter(results,
            sample,
            c('chowliu_loglik', 'prim_no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist.median[[plot.id]] = m
        plot.id = plot.id + 1
    }
    pdf(paste('plot/prim_chowliu', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()
    
    #pdf(paste('plot/prim_chowliu', branching, '_median.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist.median)
    #dev.off()
}
