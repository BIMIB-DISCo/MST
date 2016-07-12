

import.scite.output = function(dataset,
	sample.type,
	branching) {

	sample_sizes_single_cells = c(10, 25, 50, 75, 100)
	sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)

	if (sample.type == "single") {
        sample.size = sample_sizes_single_cells
    } else if (sample.type == "multiple") {
        sample.size = sample_sizes_multiple_biopses
    } else {
        stop('sample.type must be "single" or "multiple"\n')
    }

    if (! branching %in% c('low', 'medium', 'high', 'random_5', 'random_10', 'random_15', 'random_20', 'random_forest')) {
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
    			filename = paste0('scite_output/datasets/', branching,  '/', sample.type, '/', sample_size_level, '_', experiment, '_', noise, '_ml0')
                read = readLines(paste0(filename, '.gv'))
                read = gsub(' ', '', read)
                read = gsub(';', '', read)
                write(paste0(read, collapse = '\n'), file = paste0(filename, '.correct.gv'))
                readdot = read.dot(paste0(filename, '.correct.gv'))
                graph = graph.adjacency(readdot)
                matrix = get.adjacency(graph, sparse = FALSE)
                object = execution[[noise]]
                data = object$dataset
                reconstructions = object$reconstructions

                result = matrix(0, ncol = ncol(data), nrow = ncol(data))
                for (i in 1:nrow(result)) {
                    for(j in 1:ncol(result)) {
                        if (matrix[as.character(i), as.character(j)] == 1) {
                            result[i,j] = 1
                        }
                    }
                }

                true_tree = object$true_tree


                statistics = getStats(true_tree, result)
                object$reconstructions$scite$no.reg.adj$scite_no_reg = result
                object$reconstructions$scite$no.reg.res = statistics
                dataset[[sample, experiment]][[noise]] = object
    		}
        }
    }
    return(dataset)
}

