

create.scite.input = function(dataset,
    sample.type,
    branching,
    betasd,
    seed = 12345,
    pass.error.rates = TRUE) {

    if (! branching %in% c('low', 'medium', 'high', 'random_5', 'random_10', 'random_15', 'random_20', 'random_forest')) {
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
    
    count = 1
    sample.count = 1
    for (sample in 1:nrow(dataset)) {
        script = ''
        for (experiment in 1:ncol(dataset)) {
            execution = dataset[[sample, experiment]]
            for (noise in ls(execution)) {
                noise_level = as.numeric(noise)
                epos_level = execution[[noise]]$epos
                if (epos_level == 0 || !pass.error.rates) {
                    epos_level = '0.00000000000001'
                }
                eneg_level = execution[[noise]]$eneg
                if (eneg_level == 0 || !pass.error.rates) {
                    eneg_level = '0.00000000000001'
                }
                sample_size_level = rownames(dataset)[sample]
                filename = paste0('datasets/', branching,  '/', sample.type, '/', sample_size_level, '_', experiment, '_', noise, '.csv')
                genotype = t(execution[[noise]]$dataset)
                nsample = ncol(genotype)
                nmuts = nrow(genotype)
                seed_level = round(runif(1) * 10000, 0)
                write.table(genotype, file = paste0('scite_input/', filename), quote = FALSE, sep = ' ', col.names = FALSE, row.names = FALSE)
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
