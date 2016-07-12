library(ggplot2)
library(Rmisc)

source('../giulio.plot.R')
load('RData/results.values.low.RData')
results.values.low = results.values
load('RData/results.values.medium.RData')
results.values.medium = results.values
load('RData/results.values.high.RData')
results.values.high = results.values

load('RData/results.low.RData')
results.low = results
load('RData/results.medium.RData')
results.medium = results
load('RData/results.high.RData')
results.high = results

for (branching in c('low', 'medium', 'high')) {
    if (branching == 'low') {
        branching = 'lowr'
        results.values = results.values.low
        results = results.low
    } else if (branching == 'medium') {
        branching = 'mediumr'
        results.values = results.values.medium
        results = results.medium
    } else if (branching == 'high') {
        branching = 'highr'
        results.values = results.values.high
        results = results.high
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
    pdf(paste0('plot/capri_caprese_', branching, '_random_cols.pdf'), height = 14, width = 11)
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
            sample,
            c('gabow_entropy.no.reg', 'gabow_pmi.no.reg', 'gabow_cpmi.no.reg',
            'gabow_mi.no.reg', 'gabow_no_rising_no.raising.entropy.no.reg',
            'gabow_no_rising_no.raising.pmi.no.reg', 'gabow_no_rising_no.raising.cpmi.no.reg',
            'gabow_no_rising_no.raising.mi.no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist.median[[plot.id]] = m
        plot.id = plot.id + 1
    }
    pdf(paste0('plot/gabow_', branching, '_random_cols.pdf'), height = 14, width = 11)
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
        p = dotplotter(results.values, sample, c('edmonds_entropy.no.reg',
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
    pdf(paste0('plot/edmonds_', branching, '_random_cols.pdf'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()

    #pdf(paste('plot/edmonds', branching, '_median.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist.median)
    #dev.off()

    # SCITE
    #plotlist = list()
    #plotlist.median = list()
    #plot.id = 1

    #for (sample in c(10, 25, 50, 75, 100)) {
    #    p = dotplotter(results.values, sample, c('scite_no.reg'), 
    #        branching,
    #        paste('SAMPLE SIZE = ', sample))
    #    
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
            c('chowliu_loglik', 'prim_no.reg', 'scite_no.reg'),
            branching,
            paste('SAMPLE SIZE = ', sample))
        plotlist.median[[plot.id]] = m
        plot.id = plot.id + 1
    }
    pdf(paste0('plot/prim_chowliu_', branching, '_random_cols.pdf'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()
    
    #pdf(paste('plot/prim_chowliu', branching, '_median.pdf', sep='_'), height = 14, width = 11)
    #multiplot(plotlist = plotlist.median)
    #dev.off()    
}
