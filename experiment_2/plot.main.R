library(ggplot2)
library(Rmisc)

source('../giulio.plot.R')
load('RData/results.values.low.RData')
load('RData/results.values.medium.RData')
load('RData/results.values.high.RData')


for (branching in c('low', 'medium', 'high')) {
    if (branching == 'low') {
        results.values = results.values.low
    } else if (branching == 'medium') {
        results.values = results.values.medium
    } else if (branching == 'hing') {
        results.values = results.values.high
    }


    # CAPRI CAPRESE
    plotlist = list()
    plot.id = 1

    cat('capri caprese', branching, '\n')

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, 
            sample,
            c('capri_bic', 'capri_aic', 'capri_loglik', 'caprese_no.reg'),
            branching,
            paste('SAMPLE SIZE = ', sample))
        
        plotlist[[plot.id]] = p
        plot.id = plot.id + 1
    }
    pdf(paste('plot/capri_caprese', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()


    # GABOW
    plotlist = list()
    plot.id = 1

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, sample, c('gabow_entropy.no.reg', 'gabow_pmi.no.reg', 'gabow_cpmi.no.reg',
            'gabow_mi.no.reg', 'gabow_no_rising_no.raising.entropy.no.reg',
            'gabow_no_rising_no.raising.pmi.no.reg', 'gabow_no_rising_no.raising.cpmi.no.reg',
            'gabow_no_rising_no.raising.mi.no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))

        plotlist[[plot.id]] = p
        plot.id = plot.id + 1
    }
    pdf(paste('plot/gabow', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()

    # EDMONDS
    plotlist = list()
    plot.id = 1

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, sample, c('edmonds_entropy.no.reg',
            'edmonds_pmi.no.reg',
            'edmonds_cpmi.no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))

        plotlist[[plot.id]] = p
        plot.id = plot.id + 1
    }
    pdf(paste('plot/edmonds', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()

    # SCITE
    plotlist = list()
    plot.id = 1

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, sample, c('scite_no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))

        plotlist[[plot.id]] = p
        plot.id = plot.id + 1
    }
    pdf(paste('plot/scite', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()

    # PRIM CHOW LIU
    plotlist = list()
    plot.id = 1

    for (sample in c(10, 25, 50, 75, 100)) {
        p = dotplotter(results.values, sample, c('chowliu_loglik', 'prim_no.reg'), 
            branching,
            paste('SAMPLE SIZE = ', sample))

        plotlist[[plot.id]] = p
        plot.id = plot.id + 1
    }
    pdf(paste('plot/prim_chowliu', branching, '.pdf', sep='_'), height = 14, width = 11)
    multiplot(plotlist = plotlist)
    dev.off()
    
}
