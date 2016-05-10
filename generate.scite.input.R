

create.scite.input = function(dataset,
    sample.type,
    branching,
    betasd,
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

    if (! branching %in% c('low', 'medium', 'high', 'random_5', 'random_10', 'random_15', 'random_20')) {
        stop('branching must be "low", "medium" or "high"')
    } 

    if (! dir.exists('scite_input')) {
        dir.create('scite_input')
    }

    if (! dir.exists('scite_input/datasets/')) {
        dir.create('scite_input/datasets/')
    }

    if (! dir.exists(paste0('scite_input/datasets/', branching))) {
        dir.create(paste0('scite_input/datasets/', branching))
    }

    if (! dir.exists(paste0('scite_input/datasets/', branching, '/', sample.type))) {
        dir.create(paste0('scite_input/datasets/', branching, '/', sample.type))
    }

    set.seed(seed)
    runif

    count = 1
    sample.count = 1
    for (sample in 1:nrow(dataset)) {
        script = ''
        for (experiment in 1:ncol(dataset)) {
            execution = dataset[[sample, experiment]]
            for (noise in ls(execution)) {
                noise_level = as.numeric(noise)
                epos_level = epos[noise_level]
                if (epos_level == 0) {
                    epos_level = '0.00000000000001'
                }
                eneg_level = eneg[noise_level]
                eneg_level = generate.beta(eneg_level, betasd)
                sample_size_level = sample.size[sample]
                filename = paste0('scite_input/datasets/', branching,  '/', sample.type, '/', sample_size_level, '_', experiment, '_', noise, '.csv')
                genotype = t(execution[[noise]]$dataset)
                nsample = ncol(genotype)
                nmuts = nrow(genotype)
                seed_level = round(runif(1) * 10000, 0)
                write.table(genotype, file = filename, quote = FALSE, sep = ' ', col.names = FALSE, row.names = FALSE)
                script = paste0(script, 'printf "EXEC: ', count, '/800 - sample: ', sample, ' experiment: ', experiment, ' noise lev: ', noise, ' \\n" \n')
                script = paste0(script, '../scite -i ', filename, ' -n ', nmuts, ' -m ', nsample, ' -r 1 -l 100000 -fd ', epos_level, ' -ad ', eneg_level, ' -seed ', seed_level, ' -max_treelist_size 1', ' \n')
                count = count + 1
            }
        }

        script.filename = paste0('scite_input/scite.script.', sample.type, '.', branching, '.', sample.count, '.sh')
        write(script, script.filename)
        sample.count = sample.count + 1
    }
}

generate.beta <- function(mean, sd, min = 0, max = 1) {
    if (sd == 0) {
        return(mean)
    }
    beta = rnorm(1, mean, sd)
    while (beta <= min || beta >= max) {
        beta = rnorm(1, mean, sd)
    }
    return(beta)
}
