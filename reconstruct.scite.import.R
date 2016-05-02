

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

