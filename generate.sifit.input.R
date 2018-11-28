##############################################################################
###
### Generate SiFit Input
###
##############################################################################
## Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

branching.choices = c('low',
                      'medium',
                      'high',
                      'random_5',
                      'random_10',
                      'random_15',
                      'random_20',
                      'random_forest',
                      'clean',
                      'convergent',
                      'random_columns',
                      'random_forest_fixed',
                      'mini_01',
                      'missing',
                      'mini_02',
                      'mini_03',
                      'mini_04',
                      'mini_05')


create.sifit.input = function(dataset,
                              sample.type,
                              branching,
                              betasd,
                              numsample = NA,
                              seed = 12345,
                              pass.error.rates = TRUE,
                              steps = 100000) {
    
    options(scipen = 999)

    if (! branching %in% branching.choices) {
        stop('branching must be "low", "medium" or "high"')
    } 

    if (! dir.exists('sifit_input')) {
        dir.create('sifit_input')
    }

    if (! dir.exists('sifit_input/datasets/')) {
        dir.create('sifit_input/datasets/')
    }

    if (! dir.exists(paste0('sifit_input/datasets/', branching))) {
        dir.create(paste0('sifit_input/datasets/', branching))
    }

    if (! dir.exists(paste0('sifit_input/datasets/', branching, '/', sample.type))) {
        dir.create(paste0('sifit_input/datasets/', branching, '/', sample.type))
    }

    set.seed(seed)
    
    count = 1
    sample.count = 1
    for (sample in 1:nrow(dataset)) {
        script = ''
        for (experiment in 1:ncol(dataset)) {
            execution = dataset[[sample, experiment]]
            for (noise in 1:length(execution)) {
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
                if (!is.na(numsample) && numsample != sample_size_level) {
                    break;
                }
                filename = paste0('datasets/', branching,  '/', sample.type, '/',
                                  sample_size_level, '_', experiment, '_', noise,
                                  '.csv')
                genotype = t(execution[[noise]]$dataset)
                if (branching == 'missing') {
                    for (i in 1:nrow(genotype)) {
                        for(j in 1:ncol(genotype)) {
                            if (is.na(genotype[[i, j]])) {
                                genotype[[i, j]] = 3
                            }
                        }
                    }
                }
                nsample = ncol(genotype)
                nmuts = nrow(genotype)
                genotype = cbind(1:nrow(genotype), genotype)
                seed_level = round(runif(1) * 10000, 0)
                write.table(genotype,
                            file = paste0('sifit_input/', filename),
                            quote = FALSE,
                            sep = ' ',
                            col.names = FALSE,
                            row.names = FALSE)
                script = paste0(script,
                                'printf "EXEC: ', count,
                                '/800 - sample: ', sample,
                                ' experiment: ', experiment,
                                ' noise lev: ', noise,
                                ' \\n" \n')
                script = paste0(script, 'java -jar ../SiFit.jar -df 1 -ipMat ', filename,
                                ' -n ', nmuts,
                                ' -m ', nsample,
                                ' -r 1 -iter ', steps,
                                '-fp ', epos_level,
                                ' -fn ', eneg_level,
                                ' \n')
                count = count + 1
            }
        }

        script.filename = paste0('sifit_input/sifit.script.', sample.type, '.', branching, '.', sample.count, '.sh')
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


### end of file -- generate.sifit.input.R
