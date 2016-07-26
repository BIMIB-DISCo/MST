library(ggplot2)
library(Rmisc)

source('../giulio.plot.R')
e <- new.env()

for (type in c('accuracy', 'hamming_distance')) {
    for (branching in c('low', 'medium', 'high')) {
        load(paste0('RData/results.values.', type, '.', branching, '.RData'), envir = e)
        results.values = e$results.values
        load(paste0('RData/results.', type, '.', branching, '.RData'), envir = e)
        results = e$results

        # CAPRI CAPRESE
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        cat('capri caprese', branching, '\n')

        for (sample in c(5, 7, 10, 20, 50)) {
            p = dotplotter(results.values, 
                sample,
                c('capri_bic', 'capri_aic', 'capri_loglik', 'caprese_no.reg', 'scite_no.reg'),
                branching,
                type,
                paste('SAMPLE SIZE = ', sample),
                sample.type = 'multiple')
            
            plotlist[[plot.id]] = p

            m = medianplotter(results,
                sample,
                c('capri_bic', 'capri_aic', 'capri_loglik', 'caprese_no.reg'),
                branching,
                paste('SAMPLE SIZE = ', sample))
            plotlist.median[[plot.id]] = m
            plot.id = plot.id + 1
        }
        pdf(paste('plot/capri_caprese', type, branching, '.pdf', sep='_'), height = 14, width = 11)
        multiplot(plotlist = plotlist)
        dev.off()

        #pdf(paste('plot/capri_caprese', branching, '_median.pdf', sep='_'), height = 14, width = 11)
        #multiplot(plotlist = plotlist.median)
        #dev.off()

        # GABOW
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        for (sample in c(5, 7, 10, 20, 50)) {
            p = dotplotter(results.values, sample, c('gabow_entropy.no.reg', 'gabow_pmi.no.reg', 'gabow_cpmi.no.reg',
                'gabow_mi.no.reg', 'scite_no.reg'), 
                branching,
                type,
                paste('SAMPLE SIZE = ', sample),
                sample.type = 'multiple')

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
        pdf(paste('plot/gabow', type, branching, '.pdf', sep='_'), height = 14, width = 11)
        multiplot(plotlist = plotlist)
        dev.off()

        #pdf(paste('plot/gabow', branching, '_median.pdf', sep='_'), height = 14, width = 11)
        #multiplot(plotlist = plotlist.median)
        #dev.off()

        # EDMONDS
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        for (sample in c(5, 7, 10, 20, 50)) {
            p = dotplotter(results.values, sample, c('edmonds_entropy.no.reg',
                'edmonds_pmi.no.reg',
                'edmonds_cpmi.no.reg',
                'scite_no.reg'), 
                branching,
                type,
                paste('SAMPLE SIZE = ', sample),
                sample.type = 'multiple')

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
        pdf(paste('plot/edmonds', type, branching, '.pdf', sep='_'), height = 14, width = 11)
        multiplot(plotlist = plotlist)
        dev.off()

        #pdf(paste('plot/edmonds', branching, '_median.pdf', sep='_'), height = 14, width = 11)
        #multiplot(plotlist = plotlist.median)
        #dev.off()

        # PRIM CHOW LIU
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        for (sample in c(5, 7, 10, 20, 50)) {
            p = dotplotter(results.values, sample, c('chowliu_loglik', 'prim_no.reg', 'scite_no.reg'), 
                branching,
                type,
                paste('SAMPLE SIZE = ', sample),
                sample.type = 'multiple')

            plotlist[[plot.id]] = p

            m = medianplotter(results,
                sample,
                c('chowliu_loglik', 'prim_no.reg', 'scite_no.reg'),
                branching,
                paste('SAMPLE SIZE = ', sample))
            plotlist.median[[plot.id]] = m
            plot.id = plot.id + 1
        }
        pdf(paste('plot/prim_chowliu', type, branching, '.pdf', sep='_'), height = 14, width = 11)
        multiplot(plotlist = plotlist)
        dev.off()
        
        #pdf(paste('plot/prim_chowliu', branching, '_median.pdf', sep='_'), height = 14, width = 11)
        #multiplot(plotlist = plotlist.median)
        #dev.off()


        # ALL
        plotlist = list()
        plotlist.median = list()
        plot.id = 1

        for (sample in c(5, 7, 10, 20, 50)) {
            p = dotplotter(results.values, sample, c('capri_bic', 'capri_aic',
                'caprese_no.reg', 'edmonds_entropy.no.reg',
                'edmonds_pmi.no.reg', 'edmonds_cpmi.no.reg',
                'gabow_entropy.no.reg', 'gabow_pmi.no.reg', 'gabow_cpmi.no.reg',
                'gabow_mi.no.reg', 'chowliu_loglik', 'prim_no.reg',
                'scite_no.reg'), 
                branching,
                type,
                paste('SAMPLE SIZE = ', sample),
                sample.type = 'multiple')

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
        pdf(paste('plot/all', type, branching, '.pdf', sep='_'), height = 18, width = 22)
        multiplot(plotlist = plotlist)
        dev.off()



    }
}

