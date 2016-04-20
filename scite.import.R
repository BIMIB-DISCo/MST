library(igraph)
library(sna)
library(Rgraphviz)

source('statistics.R')
source('performance.plot.R')

number_experiments = 100
my_experiments = 1:number_experiments
names(my_experiments) = paste("Experiment",my_experiments)
sample_levels = 10
noise_levels = 8
my_algorithms = c("capri", "caprese", "edmonds", "chowliu", "prim", "scite")
my_regularizators = c("no.reg.res","loglik.res","aic.res","bic.res")

# set the true trees
low_true_tree = array(0,c(6,6))
low_true_tree[1,2] = 1
low_true_tree[1,3] = 1
low_true_tree[1,4] = 1
low_true_tree[2,5] = 1
low_true_tree[3,6] = 1

medium_true_tree = array(0,c(11,11))
medium_true_tree[1,2] = 1
medium_true_tree[1,3] = 1
medium_true_tree[1,4] = 1
medium_true_tree[2,5] = 1
medium_true_tree[4,6] = 1
medium_true_tree[5,7] = 1
medium_true_tree[5,8] = 1
medium_true_tree[6,9] = 1
medium_true_tree[6,10] = 1
medium_true_tree[6,11] = 1

high_true_tree = array(0,c(17,17))
high_true_tree[1,2] = 1
high_true_tree[1,3] = 1
high_true_tree[1,4] = 1
high_true_tree[2,5] = 1
high_true_tree[4,6] = 1
high_true_tree[5,7] = 1
high_true_tree[5,8] = 1
high_true_tree[6,9] = 1
high_true_tree[6,10] = 1
high_true_tree[6,11] = 1
high_true_tree[7,12] = 1
high_true_tree[9,13] = 1
high_true_tree[9,14] = 1
high_true_tree[9,15] = 1
high_true_tree[11,16] = 1
high_true_tree[11,17] = 1




import.scite.output = function(dataset,
	sample.type,
	branching,
    true_tree,
	seed = 12345) {

	e_pos_single_cells = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
	e_neg_single_cells = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
	sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
	e_pos_multiple_biopses = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
	e_neg_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
	sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)

	if (sample.type == "single") {
        epos = e_pos_single_cells
        eneg = e_neg_single_cells
        sample.size = sample_sizes_single_cells
    } else if (sample.type == "multiple") {
        epos = e_pos_multiple_biopses
        eneg = e_neg_multiple_biopses
        sample.size = sample_sizes_multiple_biopses
    } else {
        stop('sample.type must be "single" or "multiple"\n')
    }

    if (! branching %in% c('low', 'medium', 'high')) {
        stop('branching must be "low", "medium" or "high"')
    } 

    if (! dir.exists('RData')) {
    	stop('no buono')
    }

    for (sample in 1:nrow(dataset)) {
    	for (experiment in 1:ncol(dataset)) {
    		execution = dataset[[sample, experiment]]
    		for (noise in ls(execution)) {
    			noise_level = as.numeric(noise)
    			sample_size_level = sample.size[sample]
    			filename = paste0('scite_output/', branching,  '/', sample.type, '/', sample_size_level, '_', experiment, '_', noise, '_ml0')
    			#print(paste0(filename, '.gv'))
                read = readLines(paste0(filename, '.gv'))
                read = gsub(' ', '', read)
                read = gsub(';', '', read)
                write(paste0(read, collapse = '\n'), file = paste0(filename, '.correct.gv'))
                readdot = read.dot(paste0(filename, '.correct.gv'))
                graph = graph.adjacency(readdot)
                matrix = get.adjacency(graph, sparse = FALSE)
                object = execution[[noise]]
                data = object$data
                reconstructions = object$reconstructions
                #print(matrix)
                #print(data)
                #plot(graph)
                result = matrix(0, ncol = ncol(data), nrow = ncol(data))
                for (i in 1:nrow(result)) {
                    for(j in 1:ncol(result)) {
                        if (matrix[as.character(i), as.character(j)] == 1) {
                            result[i,j] = 1
                        }
                    }
                }
                #print(result)
                statistics = getStats(true_tree, result)
                #print(statistics)
                object$reconstructions$scite$no.reg.adj$scite_no_reg = result
                object$reconstructions$scite$no.reg.res = statistics
                dataset[[sample, experiment]][[noise]] = object
    		}
        }
    }
    return(dataset)
}

load('RData/experiments.single.cells.low.RData')
#load('RData/experiments.single.cells.medium.RData')
#load('RData/experiments.single.cells.high.RData')
load('RData/experiments.multiple.biopses.low.RData')
load('RData/experiments.multiple.biopses.medium.RData')
load('RData/experiments.multiple.biopses.high.RData')

#import.scite.output(experiments.single.cells.medium, 'single', 'medium')
#import.scite.output(experiments.single.cells.high, 'single', 'high')

experiments.single.cells.low.scite = import.scite.output(experiments.single.cells.low, 'single', 'low', low_true_tree)
save(experiments.single.cells.low.scite, file = 'RData/experiments.single.cells.low.scite.RData')
experiments.single.cells.low.scite.stats = get.stats(experiments.single.cells.low.scite,
    my_algorithms,
    my_regularizators,
    sample_levels,
    noise_levels,
    number_experiments)
save(experiments.single.biopses.low.scite.stats, file="RData/experiments.single.biopses.low.scite.stats.RData")
performance.plot(experiments.single.biopses.low.scite.stats, 'single', 'low')
compare.performance.plot(experiments.single.biopses.low.scite.stats, 'single', 'low')
compare.performance.plot.2d(experiments.single.biopses.low.scite.stats, 'single', 'low')

experiments.multiple.biopses.low.scite = import.scite.output(experiments.multiple.biopses.low, 'multiple', 'low', low_true_tree)
save(experiments.multiple.biopses.low.scite, file = 'RData/experiments.multiple.biopses.low.scite.RData')
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

experiments.multiple.biopses.medium.scite = import.scite.output(experiments.multiple.biopses.medium, 'multiple', 'medium', medium_true_tree)
save(experiments.multiple.biopses.medium.scite, file = 'RData/experiments.multiple.biopses.medium.scite.RData')
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

experiments.multiple.biopses.high.scite = import.scite.output(experiments.multiple.biopses.high, 'multiple', 'high', high_true_tree)
save(experiments.multiple.biopses.high.scite, file = 'RData/experiments.multiple.biopses.high.scite.RData')
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

