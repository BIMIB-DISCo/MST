source('statistics.plot.R')
source('statistics.compute.R')

number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
my_algorithms = c("capri","caprese","edmonds","chowliu","prim", "scite")
my_regularizators = c("no.reg.res","loglik.res","aic.res","bic.res")
sample_levels = 10
noise_levels = 8

load('RData/experiments.single.cells.low.scite.RData')
load('RData/experiments.single.cells.medium.scite.RData')
load('RData/experiments.single.cells.high.scite.RData')
load('RData/experiments.multiple.biopses.low.scite.RData')
load('RData/experiments.multiple.biopses.medium.scite.RData')
load('RData/experiments.multiple.biopses.high.scite.RData')



experiments.single.cells.low.scite.stats = get.stats(experiments.single.cells.low.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.single.cells.low.scite.stats, file="RData/experiments.single.cells.low.scite.stats.RData")
performance.plot(experiments.single.cells.low.scite.stats, 'single', 'low')
compare.performance.plot(experiments.single.cells.low.scite.stats, 'single', 'low')
compare.performance.plot.2d(experiments.single.cells.low.scite.stats, 'single', 'low')


experiments.single.cells.medium.scite.stats = get.stats(experiments.single.cells.medium.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.single.cells.medium.scite.stats, file="RData/experiments.single.cells.medium.scite.stats.RData")
performance.plot(experiments.single.cells.medium.scite.stats, 'single', 'medium')
compare.performance.plot(experiments.single.cells.medium.scite.stats, 'single', 'medium')
compare.performance.plot.2d(experiments.single.cells.medium.scite.stats, 'single', 'medium')


experiments.single.cells.high.scite.stats = get.stats(experiments.single.cells.high.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.single.cells.high.scite.stats, file="RData/experiments.single.cells.high.scite.stats.RData")
performance.plot(experiments.single.cells.high.scite.stats, 'single', 'high')
compare.performance.plot(experiments.single.cells.high.scite.stats, 'single', 'high')
compare.performance.plot.2d(experiments.single.cells.high.scite.stats, 'single', 'high')


experiments.multiple.biopses.low.scite.stats = get.stats(experiments.multiple.biopses.low.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.multiple.biopses.low.scite.stats, file="RData/experiments.multiple.biopses.low.scite.stats.RData")
performance.plot(experiments.multiple.biopses.low.scite.stats, 'multiple', 'low')
compare.performance.plot(experiments.multiple.biopses.low.scite.stats, 'multiple', 'low')
compare.performance.plot.2d(experiments.multiple.biopses.low.scite.stats, 'multiple', 'low')


experiments.multiple.biopses.medium.scite.stats = get.stats(experiments.multiple.biopses.medium.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.multiple.biopses.medium.scite.stats, file="RData/experiments.multiple.biopses.medium.scite.stats.RData")
performance.plot(experiments.multiple.biopses.medium.scite.stats, 'multiple', 'medium')
compare.performance.plot(experiments.multiple.biopses.medium.scite.stats, 'multiple', 'medium')
compare.performance.plot.2d(experiments.multiple.biopses.medium.scite.stats, 'multiple', 'medium')


experiments.multiple.biopses.high.scite.stats = get.stats(experiments.multiple.biopses.high.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.multiple.biopses.high.scite.stats, file="RData/experiments.multiple.biopses.high.scite.stats.RData")
performance.plot(experiments.multiple.biopses.high.scite.stats, 'multiple', 'high')
compare.performance.plot(experiments.multiple.biopses.high.scite.stats, 'multiple', 'high')
compare.performance.plot.2d(experiments.multiple.biopses.high.scite.stats, 'multiple', 'high')









performance.plot(experiments.single.cells.low.stats, 'single', 'low')
performance.plot(experiments.single.cells.medium.stats, 'single', 'medium')
performance.plot(experiments.single.cells.high.stats, 'single', 'high')
performance.plot(experiments.multiple.biopses.low.stats, 'multiple', 'low')
performance.plot(experiments.multiple.biopses.medium.stats, 'multiple', 'medium')
performance.plot(experiments.multiple.biopses.high.stats, 'multiple', 'high')

compare.performance.plot(experiments.single.cells.low.stats, 'single', 'low')
compare.performance.plot(experiments.single.cells.medium.stats, 'single', 'medium')
compare.performance.plot(experiments.single.cells.high.stats, 'single', 'high')
compare.performance.plot(experiments.multiple.biopses.low.stats, 'multiple', 'low')
compare.performance.plot(experiments.multiple.biopses.medium.stats, 'multiple', 'medium')
compare.performance.plot(experiments.multiple.biopses.high.stats, 'multiple', 'high')

compare.performance.plot.2d(experiments.single.cells.low.stats, 'single', 'low')
compare.performance.plot.2d(experiments.single.cells.medium.stats, 'single', 'medium')
compare.performance.plot.2d(experiments.single.cells.high.stats, 'single', 'high')
compare.performance.plot.2d(experiments.multiple.biopses.low.stats, 'multiple', 'low')
compare.performance.plot.2d(experiments.multiple.biopses.medium.stats, 'multiple', 'medium')
compare.performance.plot.2d(experiments.multiple.biopses.high.stats, 'multiple', 'high')
