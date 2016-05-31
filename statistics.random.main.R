source('statistics.plot.R')
source('statistics.compute.R')

number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
#my_algorithms = c("capri","caprese","edmonds","mle","chowliu","prim", "scite")
my_algorithms = c("edmonds","mle","scite", 'mltree')
#my_regularizators = c("no.reg.res","loglik.res","aic.res","bic.res")
my_regularizators = c("pmi.no.reg.res", "cpmi.no.reg.res", "entropy.no.reg.res", "no.reg.res")
sample_levels = 10
noise_levels = 8

load('RData/experiments.single.cells.low.scite.RData')
load('RData/experiments.single.cells.medium.scite.RData')
load('RData/experiments.single.cells.high.scite.RData')
load('RData/experiments.multiple.biopses.low.scite.RData')
load('RData/experiments.multiple.biopses.medium.scite.RData')
load('RData/experiments.multiple.biopses.high.scite.RData')


experiments.random.single.cells.5.nodes.scite.stats = get.stats(
    experiments.random.single.cells.5.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.single.cells.5.nodes.scite.stats, file='RData/experiments.random.single.cells.5.nodes.scite.stats.Rdata')

experiments.random.single.cells.10.nodes.scite.stats = get.stats(
    experiments.random.single.cells.10.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.single.cells.10.nodes.scite.stats, file='RData/experiments.random.single.cells.10.nodes.scite.stats.Rdata')


experiments.random.single.cells.15.nodes.scite.stats = get.stats(
    experiments.random.single.cells.15.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.single.cells.15.nodes.scite.stats, file='RData/experiments.random.single.cells.15.nodes.scite.stats.Rdata')


experiments.random.single.cells.20.nodes.scite.stats = get.stats(
    experiments.random.single.cells.20.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.single.cells.20.nodes.scite.stats, file='RData/experiments.random.single.cells.20.nodes.scite.stats.Rdata')


experiments.random.multiple.biopses.5.nodes.scite.stats = get.stats(
    experiments.random.multiple.biopses.5.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.multiple.biopses.5.nodes.scite.stats, file='RData/experiments.random.multiple.biopses.5.nodes.scite.stats.RData')

experiments.random.multiple.biopses.10.nodes.scite.stats = get.stats(
    experiments.random.multiple.biopses.10.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.multiple.biopses.10.nodes.scite.stats, file='RData/experiments.random.multiple.biopses.10.nodes.scite.stats.RData')


experiments.random.multiple.biopses.15.nodes.scite.stats = get.stats(
    experiments.random.multiple.biopses.15.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.multiple.biopses.15.nodes.scite.stats, file='RData/experiments.random.multiple.biopses.15.nodes.scite.stats.RData')


experiments.random.multiple.biopses.20.nodes.scite.stats = get.stats(
    experiments.random.multiple.biopses.20.nodes.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.random.multiple.biopses.20.nodes.scite.stats, file='RData/experiments.random.multiple.biopses.20.nodes.scite.stats.RData')



performance.plot(experiments.random.single.cells.5.nodes.scite.stats , 'single', 'random_5')
performance.plot(experiments.random.single.cells.10.nodes.scite.stats , 'single', 'random_10')
performance.plot(experiments.random.single.cells.15.nodes.scite.stats , 'single', 'random_15')
performance.plot(experiments.random.single.cells.20.nodes.scite.stats , 'single', 'random_20')

performance.plot(experiments.random.multiple.biopses.5.nodes.scite.stats , 'multiple', 'random_5')
performance.plot(experiments.random.multiple.biopses.10.nodes.scite.stats , 'multiple', 'random_10')
performance.plot(experiments.random.multiple.biopses.15.nodes.scite.stats , 'multiple', 'random_15')
performance.plot(experiments.random.multiple.biopses.20.nodes.scite.stats , 'multiple', 'random_20')

compare.performance.plot(experiments.random.single.cells.5.nodes.scite.stats , 'single', 'random_5')
compare.performance.plot(experiments.random.single.cells.10.nodes.scite.stats , 'single', 'random_10')
compare.performance.plot(experiments.random.single.cells.15.nodes.scite.stats , 'single', 'random_15')
compare.performance.plot(experiments.random.single.cells.20.nodes.scite.stats , 'single', 'random_20')

compare.performance.plot(experiments.random.multiple.biopses.5.nodes.scite.stats , 'multiple', 'random_5')
compare.performance.plot(experiments.random.multiple.biopses.10.nodes.scite.stats , 'multiple', 'random_10')
compare.performance.plot(experiments.random.multiple.biopses.15.nodes.scite.stats , 'multiple', 'random_15')
compare.performance.plot(experiments.random.multiple.biopses.20.nodes.scite.stats , 'multiple', 'random_20')

compare.performance.plot.2d(experiments.random.single.cells.5.nodes.scite.stats , 'single', 'random_5')
compare.performance.plot.2d(experiments.random.single.cells.10.nodes.scite.stats , 'single', 'random_10')
compare.performance.plot.2d(experiments.random.single.cells.15.nodes.scite.stats , 'single', 'random_15')
compare.performance.plot.2d(experiments.random.single.cells.20.nodes.scite.stats , 'single', 'random_20')

compare.performance.plot.2d(experiments.random.multiple.biopses.5.nodes.scite.stats , 'multiple', 'random_5')
compare.performance.plot.2d(experiments.random.multiple.biopses.10.nodes.scite.stats , 'multiple', 'random_10')
compare.performance.plot.2d(experiments.random.multiple.biopses.15.nodes.scite.stats , 'multiple', 'random_15')
compare.performance.plot.2d(experiments.random.multiple.biopses.20.nodes.scite.stats , 'multiple', 'random_20')






stop('tutto fatto')

# queste sono per ricostruire solo tronco

experiments.random.single.cells.5.nodes.stats = get.stats(
    result.random.single.cells.5.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)

experiments.random.single.cells.10.nodes.stats = get.stats(
    result.random.single.cells.10.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)

experiments.random.single.cells.15.nodes.stats = get.stats(
    result.random.single.cells.15.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)

experiments.random.single.cells.20.nodes.stats = get.stats(
    result.random.single.cells.20.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)

experiments.random.multiple.biopses.5.nodes.stats = get.stats(
    result.random.multiple.biopses.5.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)

experiments.random.multiple.biopses.10.nodes.stats = get.stats(
    result.random.multiple.biopses.10.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)

experiments.random.multiple.biopses.15.nodes.stats = get.stats(
    result.random.multiple.biopses.15.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)

experiments.random.multiple.biopses.20.nodes.stats = get.stats(
    result.random.multiple.biopses.20.nodes,
    my_algorithms, my_regularizators,
    sample_levels, noise_levels, number_experiments)





performance.plot(experiments.random.single.cells.5.nodes.stats , 'single', 'random_5')
performance.plot(experiments.random.single.cells.10.nodes.stats , 'single', 'random_10')
performance.plot(experiments.random.single.cells.15.nodes.stats , 'single', 'random_15')
performance.plot(experiments.random.single.cells.20.nodes.stats , 'single', 'random_20')

performance.plot(experiments.random.multiple.biopses.5.nodes.stats , 'multiple', 'random_5')
performance.plot(experiments.random.multiple.biopses.10.nodes.stats , 'multiple', 'random_10')
performance.plot(experiments.random.multiple.biopses.15.nodes.stats , 'multiple', 'random_15')
performance.plot(experiments.random.multiple.biopses.20.nodes.stats , 'multiple', 'random_20')

compare.performance.plot(experiments.random.single.cells.5.nodes.stats , 'single', 'random_5')
compare.performance.plot(experiments.random.single.cells.10.nodes.stats , 'single', 'random_10')
compare.performance.plot(experiments.random.single.cells.15.nodes.stats , 'single', 'random_15')
compare.performance.plot(experiments.random.single.cells.20.nodes.stats , 'single', 'random_20')

compare.performance.plot(experiments.random.multiple.biopses.5.nodes.stats , 'multiple', 'random_5')
compare.performance.plot(experiments.random.multiple.biopses.10.nodes.stats , 'multiple', 'random_10')
compare.performance.plot(experiments.random.multiple.biopses.15.nodes.stats , 'multiple', 'random_15')
compare.performance.plot(experiments.random.multiple.biopses.20.nodes.stats , 'multiple', 'random_20')

compare.performance.plot.2d(experiments.random.single.cells.5.nodes.stats , 'single', 'random_5')
compare.performance.plot.2d(experiments.random.single.cells.10.nodes.stats , 'single', 'random_10')
compare.performance.plot.2d(experiments.random.single.cells.15.nodes.stats , 'single', 'random_15')
compare.performance.plot.2d(experiments.random.single.cells.20.nodes.stats , 'single', 'random_20')

compare.performance.plot.2d(experiments.random.multiple.biopses.5.nodes.stats , 'multiple', 'random_5')
compare.performance.plot.2d(experiments.random.multiple.biopses.10.nodes.stats , 'multiple', 'random_10')
compare.performance.plot.2d(experiments.random.multiple.biopses.15.nodes.stats , 'multiple', 'random_15')
compare.performance.plot.2d(experiments.random.multiple.biopses.20.nodes.stats , 'multiple', 'random_20')
