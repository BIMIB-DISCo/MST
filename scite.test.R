

create.scite.input = function(dataset,
	sample.type,
	branching,
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

    if (! dir.exists('scite_input')) {
    	dir.create('scite_input')
    }

    if (! dir.exists(paste0('scite_input/', branching))) {
        dir.create(paste0('scite_input/', branching))
    }

    if (! dir.exists(paste0('scite_input/', branching, '/', sample.type))) {
        dir.create(paste0('scite_input/', branching, '/', sample.type))
    }

    set.seed(seed)
    runif

    script = ''
    count = 1

    for (sample in 1:nrow(dataset)) {
    	for (experiment in 1:ncol(dataset)) {
    		execution = dataset[[sample, experiment]]
    		for (noise in ls(execution)) {
    			noise_level = as.numeric(noise)
    			epos_level = epos[noise_level]
    			if (epos_level == 0) {
    				epos_level = '0.00000000000001'
    			}
    			eneg_level = eneg[noise_level]
    			if (eneg_level == 0) {
    				eneg_level = '0.00000000000001'
    			}
    			sample_size_level = sample.size[sample]
    			filename = paste0('scite_input/', branching,  '/', sample.type, '/', sample_size_level, '_', experiment, '_', noise, '.csv')
    			genotype = t(execution[[noise]]$dataset)
    			nsample = ncol(genotype)
    			nmuts = nrow(genotype)
    			seed_level = round(runif(1) * 10000, 0)
    			write.table(genotype, file = filename, quote = FALSE, sep = ' ', col.names = FALSE, row.names = FALSE)
    			script = paste0(script, 'printf "EXEC: ', count, '/8000 - sample: ', sample, ' experiment: ', experiment, ' noise lev: ', noise, ' \\n" \n')
    			script = paste0(script, './scite -i ', filename, ' -n ', nmuts, ' -m ', nsample, ' -r 1 -l 100000 -fd ', epos_level, ' -ad ', eneg_level, ' -seed ', seed_level, ' -max_treelist_size 1', ' \n')
    			count = count + 1
    		}
    	}
    }

    script.filename = paste0('./scite.script.', sample.type, '.', branching, '.sh')
    write(script, script.filename)
	
}

load('RData/experiments.single.cells.low.RData')
load('RData/experiments.single.cells.medium.RData')
load('RData/experiments.single.cells.high.RData')
load('RData/experiments.multiple.biopses.low.RData')
load('RData/experiments.multiple.biopses.medium.RData')
load('RData/experiments.multiple.biopses.high.RData')

if (dir.exists('scite_input')) {
	unlink('scite_input', recursive = TRUE)
	unlink('scite.script.*')
}

create.scite.input(experiments.single.cells.low, 'single', 'low')
create.scite.input(experiments.single.cells.medium, 'single', 'medium')
create.scite.input(experiments.single.cells.high, 'single', 'high')
create.scite.input(experiments.multiple.biopses.low, 'multiple', 'low')
create.scite.input(experiments.multiple.biopses.medium, 'multiple', 'medium')
create.scite.input(experiments.multiple.biopses.high, 'multiple', 'high')

