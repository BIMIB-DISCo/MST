load('RData/complete.results.RData')

library(ggplot2)


xlabels = c('EXP2' = 'Original dataset',
    'MD_1' = 'Missing data 10%',                           
    'MD_2' = 'Missing data 20%',
    'MD_3' = 'Missing data 30%',
    'MD_4' = 'Missing data 40%')

#experiment.palette = c('EXP2' = 'red',
#    'MD_1' = '#fa9fb5',                           
#    'MD_2' = '#f768a1',
#    'MD_3' = '#dd3497',
#    'MD_4' = '#ae017e')


experiment.palette = c('capri' = '#fc9272',
        'caprese' = '#67001f',                           
        'edmonds' = '#f16913',
        'gabow' = '#6baed6',
        'chowliu' = '#006d2c',
        'prim' = '#66c2a4',
        'scite' = '#ffff32')

experiment.names = c(
    'capri' = 'Capri (BIC)',
    'caprese' = 'Caprese',                           
    'edmonds' = 'Edmonds (pmi)',
    'gabow' = 'Gabow (pmi)',
    'chowliu' = 'Chow-Liu',
    'prim' = 'PRIM',
    'scite' = 'SCITE')


description = c('sensitivity' = 'Sensitivity',
    'specificity' = 'Specificity')

library(ggplot2)

p = ggplot(complete.results, aes(x = source, y = sensitivity, fill = algorithm)) +
    scale_fill_manual(values = experiment.palette, labels = experiment.names) +
    scale_x_discrete(label=xlabels) +
    ylab(paste0('Sensitivity (n = 11)')) +
    geom_boxplot(outlier.size = 0) +
    xlab('epos = 0.005; eneg = 0.05')

dev.new(width = 16, height = 5)
print(p)
dev.copy2pdf(file = 'plot/sensitivity.pdf')
dev.off()


p = ggplot(complete.results, aes(x = source, y = specificity, fill = algorithm)) +
    scale_fill_manual(values = experiment.palette, labels = experiment.names) +
    scale_x_discrete(label=xlabels) +
    ylab(paste0('Specificity (n = 11)')) +
    geom_boxplot(outlier.size = 0) +
    xlab('epos = 0.005; eneg = 0.05')

dev.new(width = 16, height = 5)
print(p)
dev.copy2pdf(file = 'plot/specificity.pdf')
dev.off()
