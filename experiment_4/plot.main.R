

### statistics

load('RData/experiments.single.cells.mini.01.scite.stats.RData')
load('RData/experiments.single.cells.mini.02.scite.stats.RData')
load('RData/experiments.single.cells.mini.03.scite.stats.RData')
load('RData/experiments.single.cells.mini.04.scite.stats.RData')
load('RData/experiments.single.cells.mini.05.scite.stats.RData')

epos = NULL
eneg = NULL

epos$mini_01 = c(0.003, 0.004, 0.005, 0.006, 0.007)
epos$mini_02 = c(0.003, 0.004, 0.005, 0.006, 0.007)
epos$mini_03 = c(0.003, 0.004, 0.005, 0.006, 0.007)
epos$mini_04 = c(0.003, 0.004, 0.005, 0.006, 0.007)
epos$mini_05 = c(0.003, 0.004, 0.005, 0.006, 0.007)

eneg$mini_01 = rep(0.03, 5) 
eneg$mini_02 = rep(0.04, 5) 
eneg$mini_03 = rep(0.05, 5) 
eneg$mini_04 = rep(0.06, 5) 
eneg$mini_05 = rep(0.07, 5) 

real.epos = 0.005
real.eneg = 0.05

real.results = NULL
load('RData/results.specificity.medium.RData')
results = results[which(results$sample == 50), ]
real.results[['specificity']] = results[which(results$noise == 2), ]

load('RData/results.sensitivity.medium.RData')
results = results[which(results$sample == 50), ]
real.results[['sensitivity']] = results[which(results$noise == 2), ]


library(ggplot2)
library(Rmisc)

source('../giulio.plot.R')
samples = c(50)
e = new.env()

for (type in c('sensitivity', 'specificity')) {
    total.results = NULL
    for (branching in c('mini_01', 'mini_02', 'mini_03', 'mini_04', 'mini_05')) {
        load(paste0('RData/results.values.', type, '.', branching, '.RData'), envir = e)
        results.values = e$results.values
        load(paste0('RData/results.', type, '.', branching, '.RData'), envir = e)
        results = e$results

        # ALL
        plotlist = list()
        plotlist.jitter = list()
        plot.id = 1

#        for (sample in samples) {
#            p = dotplotter(results.values, sample, c('capri_bic',
#                'caprese_no.reg',
#                'edmonds_pmi.no.reg',
#                'gabow_pmi.no.reg',
#                'chowliu_loglik', 'prim_no.reg',
#                'scite_no.reg'), 
#                branching,
#                type,
#                paste('SAMPLE SIZE = ', sample),
#                external.epos = epos[[branching]],
#                external.eneg = eneg[[branching]])
#
#            plotlist[[plot.id]] = p
#            
#            m = jitterplotter(results.values, sample, c('capri_bic',
#                'caprese_no.reg',
#                'edmonds_pmi.no.reg',
#                'gabow_pmi.no.reg',
#                'chowliu_loglik', 'prim_no.reg',
#                'scite_no.reg'), 
#                branching,
#                type,
#                paste('SAMPLE SIZE = ', sample),
#                external.epos = epos[[branching]],
#                external.eneg = eneg[[branching]])
#            plotlist.jitter[[plot.id]] = m
#            plot.id = plot.id + 1
#        }
#        pdf(paste('plot_reduced/mt', type, branching, '.pdf', sep='_'), height = 5, width = 15)
#        multiplot(plotlist = plotlist)
#        dev.off()
#
#        pdf(paste('plot_reduced/mt_jitter', type, branching, '.pdf', sep='_'), height = 5, width = 15)
#        multiplot(plotlist = plotlist.jitter)
#        dev.off()
#
        this.epos = epos[[branching]]
        this.eneg = eneg[[branching]]
        results$epos = this.epos
        results$eneg = this.eneg
        total.results = rbind(results, total.results)
    }
    for (algo in c('capri_bic', 'caprese_no.reg', 'edmonds_pmi.no.reg', 'gabow_pmi.no.reg',
    'chowliu_loglik', 'prim_no.reg', 'scite_no.reg')) {
        cat('\n\n')
        real.mean = real.results[[type]][real.results[[type]]$algorithm == algo,]$mean
        real.sd = real.results[[type]][real.results[[type]]$algorithm == algo,]$sd
        print(real.mean)
        print(real.sd)
        pals = c("red", "white", "red")
        breaks = c(-Inf,real.mean-(real.sd/2), real.mean+(real.sd/2),Inf)


        algo.results = total.results[which(total.results$algorithm == algo), ]
        algo.results$meanplussd = algo.results$mean + (algo.results$sd / 2)
        algo.results$block = cut(algo.results$meanplussd, breaks = breaks)
        algo.results$kldivergence = log(algo.results$sd / real.sd) + (((algo.results$sd)^2 + (real.mean - algo.results$mean)^2) / (2 * (algo.results$sd)^2)) - 0.5

        dev.new()
        p = ggplot(algo.results, aes(x = epos, y = eneg)) +
            geom_tile(aes(fill = kldivergence), colour = "white") +
            #scale_fill_viridis(name="Mean") +
            #scale_fill_gradient2(low = "blue", mid = "white", high ="red", 
            #           midpoint = real.mean, space = "Lab", guide = "colourbar") +
            #scale_fill_manual(values = pals) +
            scale_x_continuous(expand=c(0,0)) + 
            scale_y_continuous(expand=c(0,0)) +
            ggtitle(paste(algo, 'M:', round(real.mean,3)))
        print(p)
        dev.copy2pdf(file = paste0('plot_reduced/', algo, '_', type, '.pdf'))
        dev.off()
        caption = paste(type, algo, 'M:', round(real.mean,3), 'SD: ', round(real.sd, 3))
        print(caption)
        print(xtable(algo.results[,c('epos','eneg','mean','sd')], 
            digits = 3, 
            caption = caption, 
            label = paste0('table:exp4-',type, '-', algo)),
            include.rownames=FALSE)

    }

}
