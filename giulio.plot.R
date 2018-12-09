##############################################################################
###
### MST
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

giulio.plot <- function(dataset,
                        sample.type,
                        branching,
                        sample_size = NULL,
                        statistics = NULL) {

    if (is.null(sample_size)) {
        sample_sizes_single_cells = c(10, 25, 50, 75, 100)
        sample_sizes_multiple_biopses = c(5, 7, 10, 20, 50)
    } else {
        sample_sizes_single_cells = sample_size
        sample_sizes_multiple_biopses = sample_size
    }

    ##if (! branching %in% c('low', 'medium', 'high', 'random_5', 'random_10', 'random_15', 'random_20',
    ##    'lowr', 'mediumr', 'highr', 'random_forest', 'clean', 'convergent', 'random_columns', 'mini')) {
    ##    stop('branching must be "low", "medium" or "high"')
    ##} 

    ##sample_sizes_single_cells = c(10, 25, 50, 75, 100)
    ##sample_sizes_multiple_biopses = c(5, 7, 10, 20, 50)

    if (sample.type == "single") {
        sample = sample_sizes_single_cells
        ##select.sample = c(10, 25, 50, 75, 100)
    } else if (sample.type == "multiple") {
        sample = sample_sizes_multiple_biopses
        ##select.sample = c(5, 7, 10, 20, 50)
    } else {
        stop('sample.type must be "single" or "multiple"\n')
    }

    if (branching %in% c('mini_01', 'mini_02', 'mini_03', 'mini_04', 'mini_05')) {
        sample = 50
        if (sample.type == 'multiple') {
            sample = 20
        }
    }

    select.sample = sample

    if (branching == 'low') {
        nodes = 6
    } else if (branching %in% c('medium', 'clean', 'convergent', 'random_columns', 
                                'mini_01', 'mini_02', 'mini_03', 'mini_04', 'mini_05')) {
        nodes = 11
    } else if (branching == 'high') {
        nodes = 17
    } else if (branching == 'random_5') {
        nodes = 5
    } else if (branching == 'random_10') {
        nodes = 10
    } else if (branching == 'random_15') {
        nodes = 15
    } else if (branching == 'random_20' || branching == 'random_forest') {
        nodes = 20
    }


    connections = (nodes * nodes) - nodes

    if(is.null(statistics)) {
        statistics = c('hamming_distance', 'accuracy', 'sensitivity', 'specificity')
    }

    for (type in statistics) {
        cat(type, ' ')

        results = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.values = data.frame(x = NULL, stringsAsFactors = FALSE)
        ordered.regs = NULL
        
        algorithms = get(type, dataset)
        for (algorithm in names(algorithms)) {
            cat('\n  ', algorithm, ' ')
            regs = get(algorithm, algorithms)
            for (reg in names(regs)) {

                stats = get(reg, regs)

                m = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                median = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                sd = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                values = matrix(list(), nrow = nrow(stats), ncol = ncol(stats))

                for (i in 1:nrow(stats)) {
                    for (j in 1:ncol(stats)) {
                        m[i,j] = stats[[i,j]]$mean
                        median[i,j] = stats[[i,j]]$median
                        sd[i,j] = stats[[i,j]]$sd
                        values[i,j] = list(unlist(stats[[i,j]]$values))
                    }
                }

                cat('\n    ', reg)
                for (i in 1:nrow(m)) {
                    for(j in 1:ncol(m)) {
                        if (sample[j] %in% select.sample) {
                            results = rbind(results, c(i,
                                                       m[i,j],
                                                       sd[i,j],
                                                       median[i,j],
                                                       sample[j],
                                                       paste(algorithm, reg, sep = '_')), stringsAsFactors = FALSE)

                            for (value in values[[i,j]]) {
                                results.values = rbind(results.values, c(i, 
                                                                         value, 
                                                                         value / connections,
                                                                         sample[j],
                                                                         paste(algorithm, reg, sep = '_')), stringsAsFactors = FALSE)
                            }
                        }
                    }
                }
                ordered.regs = c(ordered.regs, paste(algorithm, reg))              

            }
        }

        cat('\n')

        colnames(results) = c('noise', 'mean', 'sd', 'median', 'sample', 'algorithm')
        colnames(results.values) = c('noise', 'value', 'perc', 'sample', 'algorithm')

        results$noise = as.factor(results$noise)
        results$mean = as.numeric(results$mean)
        results$sd = as.numeric(results$sd)
        results$median = as.numeric(results$median)
        results$sample = as.numeric(results$sample)
        results$algorithm = as.factor(results$algorithm)

        results.values$noise = as.factor(results.values$noise)
        results.values$value = as.numeric(results.values$value)
        results.values$perc = as.numeric(results.values$perc)
        results.values$sample = as.numeric(results.values$sample)
        results.values$algorithm = as.factor(results.values$algorithm)

        save(results, file = paste0('RData/results.', type, '.', branching, '.RData'))
        save(results.values, file = paste0('RData/results.values.', type, '.', branching, '.RData'))

    }
}




add.alpha <- function(col, alpha = 1) {
    if(missing(col))
        stop("Please provide a vector of colours.")
    r = apply(sapply(col, col2rgb)/255, 2, 
              function(x) 
                  rgb(x[1], x[2], x[3], alpha = alpha))
    return(as.vector(r))  
}




dotplotter <- function(res.values, 
                       sample,
                       algorithm,
                       branching,
                       type,
                       title,
                       sample.type = 'single',
                       external.epos = NA,
                       external.eneg = NA,
                       noise = NULL) {

    res.values = res.values[which(res.values$sample == sample), ]
    res.values = res.values[which(res.values$algorithm %in% algorithm), ]
    if (!is.null(noise)) {
        res.values = res.values[which(res.values$noise %in% noise), ]
    }

    experiment.palette = c('capri_loglik' = '##ef3b2c', 
                           'capri_aic' = '##fb6a4a',
                           'capri_bic' = '##fc9272',
                           'caprese_no.reg' = '##67001f',                           
                           'edmonds_entropy.no.reg' = '##fd8d3c',
                           'edmonds_pmi.no.reg' = '##f16913',
                           'edmonds_cpmi.no.reg' = '##d94801',
                           'gabow_entropy.no.reg' = '##9ecae1',
                           'gabow_pmi.no.reg' = '##6baed6',
                           'gabow_cpmi.no.reg' = '##4292c6',
                           'gabow_mi.no.reg' = '##2171b5',
                           'gabow_no_rising_no.raising.entropy.no.reg' = '##fa9fb5',
                           'gabow_no_rising_no.raising.pmi.no.reg' = '##f768a1',
                           'gabow_no_rising_no.raising.cpmi.no.reg' = '##dd3497',
                           'gabow_no_rising_no.raising.mi.no.reg' = '##ae017e',
                           'chowliu_loglik' = '##006d2c',
                           'prim_no.reg' = '##66c2a4',
                           'scite_no.reg' = '##ffff32')

    experiment.names = c('capri_loglik' = 'Capri (logLik)', 
                         'capri_aic' = 'Capri (AIC)',
                         'capri_bic' = 'Capri (BIC)',
                         'caprese_no.reg' = 'Caprese',                           
                         'edmonds_entropy.no.reg' = 'Edmonds (entropy)',
                         'edmonds_pmi.no.reg' = 'Edmonds (pmi)',
                         'edmonds_cpmi.no.reg' = 'Edmonds (cpmi)',
                         'gabow_entropy.no.reg' = 'Gabow (entropy)',
                         'gabow_pmi.no.reg' = 'Gabow (pmi)',
                         'gabow_cpmi.no.reg' = 'Gabow (cpmi)',
                         'gabow_mi.no.reg' = 'Gabow (mi)',
                         'gabow_no_rising_no.raising.entropy.no.reg' = 'Gabow (entropy no PR)',
                         'gabow_no_rising_no.raising.pmi.no.reg' = 'Gabow (pmi no PR)',
                         'gabow_no_rising_no.raising.cpmi.no.reg' = 'Gabow (cpmi no PR)',
                         'gabow_no_rising_no.raising.mi.no.reg' = 'Gabow (mi no PR)',
                         'chowliu_loglik' = 'Chow - Liu',
                         'prim_no.reg' = 'PRIM',
                         'scite_no.reg' = 'SCITE')

    if (branching == 'low') {
        nodes = 6
        cut = 3
    } else if (branching %in% c('medium', 'clean', 'mini_01', 'mini_02',
                                'mini_03', 'mini_04', 'mini_05')) {
        nodes = 11
        cut = 4
    } else if (branching == 'convergent') {
        nodes = 8
        cut = 4
    } else if (branching == 'high') {
        nodes = 17
        cut = 8
    } else if (branching == 'random_5') {
        nodes = 5
        cut = 3
    } else if (branching == 'random_10') {
        nodes = 10
        cut = 5
    } else if (branching == 'random_15') {
        nodes = 15
        cut = 7
    } else if (branching %in% c('random_20', 'large_20')) {
        nodes = 20
        cut = 11
    } else if (branching == 'lowr') {
        nodes = 7
        cut = 3
    } else if (branching %in% c('mediumr', 'random_columns')) {
        nodes = 13
        cut = 4
    } else if (branching == 'highr') {
        nodes = 21
        cut = 8
    } else if (branching == 'random_forest') {
        nodes = 20
        cut = 11
    }

    if (type == 'hamming_distance') {
        description = 'Hamming Distance'
    } else if (type == 'accuracy') {
        description = 'Accuracy'
    } else if (type == 'sensitivity') {
        description = 'Sensitivity'
    } else if (type == 'specificity') {
        description = 'Specificity'
    } else if (type == 'elasped_time') {
        description = 'Time'
    }

    cat(title, '\n')

    epos = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
    eneg = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)

    if (sample.type == 'multiple') {
        epos = c(0.000, 0.05, 0.1, 0.15, 0.20)
        eneg = epos
    }
    if (branching %in% c('mini_01', 'mini_02',
                         'mini_03', 'mini_04', 'mini_05', 'large_20')) {
        epos = external.epos
        eneg = external.eneg
        if (sample.type == 'multiple') {
            epos = eneg
        }
    }

    if (!is.null(noise)) {
        epos = epos[noise]
        eneg = eneg[noise]
    }

    xlabels = paste(format(epos, scientific = TRUE), format(eneg, scientific = TRUE), sep = ',')
    connections = nodes * (nodes - 1)

    p = ggplot(res.values, aes(x = noise, y = value, fill = algorithm)) + 
        ##p = ggplot(res.values, aes(x = reorder(noise, value, FUN = median), y = value, fill = algorithm)) + 
        scale_fill_manual(values = experiment.palette, labels = experiment.names) +
        scale_x_discrete(label = xlabels) +
        ##expand_limits(y = c(0,connections/cut)) +
        ##scale_y_continuous(label = function(x){paste0(x, ' (', as.integer(x/connections * 100), '%)')}) +
        ylab(paste0(description, ' (n = ', nodes, ')')) +
        xlab('False Positives, False Negatives (rates)') +
        geom_boxplot(outlier.size = 0) +
        ##geom_jitter(pch = 21, position = position_jitterdodge(jitter.height = 0.2, jitter.width = 2)) +
        ggtitle(title)
    
    return(p)
}



## median

medianplotter <- function(res, sample, algorithm, branching, title) {

    res = res[which(res$sample == sample), ]
    res = res[which(res$algorithm %in% algorithm), ]

    experiment.palette = c('capri_loglik' = '#ef3b2c', 
                           'capri_aic' = '#fb6a4a',
                           'capri_bic' = '#fc9272',
                           'caprese_no.reg' = '#67001f',                           
                           'edmonds_entropy.no.reg' = '#fd8d3c',
                           'edmonds_pmi.no.reg' = '#f16913',
                           'edmonds_cpmi.no.reg' = '#d94801',
                           'gabow_entropy.no.reg' = '#9ecae1',
                           'gabow_pmi.no.reg' = '#6baed6',
                           'gabow_cpmi.no.reg' = '#4292c6',
                           'gabow_mi.no.reg' = '#2171b5',
                           'gabow_no_rising_no.raising.entropy.no.reg' = '#fa9fb5',
                           'gabow_no_rising_no.raising.pmi.no.reg' = '#f768a1',
                           'gabow_no_rising_no.raising.cpmi.no.reg' = '#dd3497',
                           'gabow_no_rising_no.raising.mi.no.reg' = '#ae017e',
                           'chowliu_loglik' = '#006d2c',
                           'prim_no.reg' = '#66c2a4',
                           'scite_no.reg' = '#ffff32')

    experiment.names = c('capri_loglik' = 'Capri (logLik)', 
                         'capri_aic' = 'Capri (AIC)',
                         'capri_bic' = 'Capri (BIC)',
                         'caprese_no.reg' = 'caprese',                           
                         'edmonds_entropy.no.reg' = 'entropy',
                         'edmonds_pmi.no.reg' = 'PMI',
                         'edmonds_cpmi.no.reg' = 'CPMI',
                         'gabow_entropy.no.reg' = 'entropy',
                         'gabow_pmi.no.reg' = 'PMI',
                         'gabow_cpmi.no.reg' = 'CPMI',
                         'gabow_mi.no.reg' = 'MI',
                         'gabow_no_rising_no.raising.entropy.no.reg' = 'entropy NR',
                         'gabow_no_rising_no.raising.pmi.no.reg' = 'PMI NR',
                         'gabow_no_rising_no.raising.cpmi.no.reg' = 'CPMI NR',
                         'gabow_no_rising_no.raising.mi.no.reg' = 'MI NR',
                         'chowliu_loglik' = 'chowliu',
                         'prim_no.reg' = 'prim',
                         'scite_no.reg' = 'SCITE')

    if (branching == 'low') {
        nodes = 6
        cut = 3
    } else if (branching %in% c('medium', 'clean', 'mini_01', 'mini_02',
                                'mini_03', 'mini_04', 'mini_05')) {
        nodes = 11
        cut = 4
    } else if (branching == 'convergent') {
        nodes = 8
        cut = 4
    } else if (branching == 'high') {
        nodes = 17
        cut = 8
    } else if (branching == 'random_5') {
        nodes = 5
        cut = 3
    } else if (branching == 'random_10') {
        nodes = 10
        cut = 5
    } else if (branching == 'random_15') {
        nodes = 15
        cut = 7
    } else if (branching == 'random_20') {
        nodes = 20
        cut = 11
    } else if (branching == 'lowr') {
        nodes = 7
        cut = 3
    } else if (branching %in% c('mediumr', 'random_columns')) {
        nodes = 13
        cut = 4
    } else if (branching == 'highr') {
        nodes = 21
        cut = 8
    } else if (branching == 'random_forest') {
        nodes = 20
        cut = 11
    }

    cat(title, '\n')

    epos = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
    eneg = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
    xlabels = paste(format(epos, scientific = TRUE), format(eneg, scientific = TRUE), sep = ',')
    connections = nodes * (nodes - 1)

    p = NULL
    p = ggplot(res, aes(x = noise, y = median, group = algorithm, colour = algorithm)) +
        scale_color_manual(values = experiment.palette, labels = experiment.names) +
        geom_line() +
        geom_point() +
        scale_x_discrete(label = xlabels) +
        expand_limits(y = c(0,nodes/2)) +
        ##scale_y_continuous(label = function(x){paste0(round(x/connections * 100, 2), '%')}) +
        ylab(paste0('Hamming distance (n = ', nodes, ')')) +
        xlab('False Positives, False Negatives (rates)') +
        ##geom_jitter(pch = 21, position = position_jitterdodge(jitter.height = 0.2, jitter.width = 2))
        ggtitle(title)
    
    return(p)
}


jitterplotter <- function(res.values,
                          sample,
                          algorithm,
                          branching,
                          type, title,
                          sample.type = 'single',
                          external.epos = NA,
                          external.eneg = NA,
                          noise = NULL) {

    res.values = res.values[which(res.values$sample == sample), ]
    res.values = res.values[which(res.values$algorithm %in% algorithm), ]
    if (!is.null(noise)) {
        res.values = res.values[which(res.values$noise %in% noise), ]
    }

    experiment.palette = c('capri_loglik' = '#ef3b2c', 
                           'capri_aic' = '#fb6a4a',
                           'capri_bic' = '#fc9272',
                           'caprese_no.reg' = '#67001f',                           
                           'edmonds_entropy.no.reg' = '#fd8d3c',
                           'edmonds_pmi.no.reg' = '#f16913',
                           'edmonds_cpmi.no.reg' = '#d94801',
                           'gabow_entropy.no.reg' = '#9ecae1',
                           'gabow_pmi.no.reg' = '#6baed6',
                           'gabow_cpmi.no.reg' = '#4292c6',
                           'gabow_mi.no.reg' = '#2171b5',
                           'gabow_no_rising_no.raising.entropy.no.reg' = '#fa9fb5',
                           'gabow_no_rising_no.raising.pmi.no.reg' = '#f768a1',
                           'gabow_no_rising_no.raising.cpmi.no.reg' = '#dd3497',
                           'gabow_no_rising_no.raising.mi.no.reg' = '#ae017e',
                           'chowliu_loglik' = '#006d2c',
                           'prim_no.reg' = '#66c2a4',
                           'scite_no.reg' = '#ffff32')

    experiment.names = c('capri_loglik' = 'Capri (logLik)', 
                         'capri_aic' = 'Capri (AIC)',
                         'capri_bic' = 'Capri (BIC)',
                         'caprese_no.reg' = 'Caprese',                           
                         'edmonds_entropy.no.reg' = 'Edmonds (entropy)',
                         'edmonds_pmi.no.reg' = 'Edmonds (pmi)',
                         'edmonds_cpmi.no.reg' = 'Edmonds (cpmi)',
                         'gabow_entropy.no.reg' = 'Gabow (entropy)',
                         'gabow_pmi.no.reg' = 'Gabow (pmi)',
                         'gabow_cpmi.no.reg' = 'Gabow (cpmi)',
                         'gabow_mi.no.reg' = 'Gabow (mi)',
                         'gabow_no_rising_no.raising.entropy.no.reg' = 'Gabow (entropy no PR)',
                         'gabow_no_rising_no.raising.pmi.no.reg' = 'Gabow (pmi no PR)',
                         'gabow_no_rising_no.raising.cpmi.no.reg' = 'Gabow (cpmi no PR)',
                         'gabow_no_rising_no.raising.mi.no.reg' = 'Gabow (mi no PR)',
                         'chowliu_loglik' = 'Chow - Liu',
                         'prim_no.reg' = 'PRIM',
                         'scite_no.reg' = 'SCITE')

    if (branching == 'low') {
        nodes = 6
        cut = 3
    } else if (branching %in% c('medium', 'clean', 'mini_01', 'mini_02',
                                'mini_03', 'mini_04', 'mini_05')) {
        nodes = 11
        cut = 4
    } else if (branching == 'convergent') {
        nodes = 8
        cut = 4
    } else if (branching == 'high') {
        nodes = 17
        cut = 8
    } else if (branching == 'random_5') {
        nodes = 5
        cut = 3
    } else if (branching == 'random_10') {
        nodes = 10
        cut = 5
    } else if (branching == 'random_15') {
        nodes = 15
        cut = 7
    } else if (branching == 'random_20') {
        nodes = 20
        cut = 11
    } else if (branching == 'lowr') {
        nodes = 7
        cut = 3
    } else if (branching %in% c('mediumr', 'random_columns')) {
        nodes = 13
        cut = 4
    } else if (branching == 'highr') {
        nodes = 21
        cut = 8
    } else if (branching == 'random_forest') {
        nodes = 20
        cut = 11
    }

    if (type == 'hamming_distance') {
        description = 'Hamming Distance'
    } else if (type == 'accuracy') {
        description = 'Accuracy'
    } else if (type == 'sensitivity') {
        description = 'Sensitivity'
    } else if (type == 'specificity') {
        description = 'Specificity'
    } else if (type == 'elasped_time') {
        description = 'Time'
    }

    cat(title, '\n')

    epos = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
    eneg = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)

    if (sample.type == 'multiple') {
        epos = c(0.000, 0.05, 0.1, 0.15, 0.20)
        eneg = epos
    }

    if (branching %in% c('mini_01', 'mini_02',
                         'mini_03', 'mini_04', 'mini_05')) {
        epos = external.epos
        eneg = external.eneg
        if (sample.type == 'multiple') {
            epos = eneg
        }
    }


    if (!is.null(noise)) {
        epos = epos[noise]
        eneg = eneg[noise]
    }

    xlabels = paste(format(epos, scientific = TRUE), format(eneg, scientific = TRUE), sep = ',')
    connections = nodes * (nodes - 1)

    p = ggplot(res.values, aes(x = noise, y = value, fill = algorithm)) + 
        ##p = ggplot(res.values, aes(x = reorder(noise, value, FUN = median), y = value, fill = algorithm)) + 
        scale_fill_manual(values = experiment.palette, labels = experiment.names) +
        scale_x_discrete(label = xlabels) +
        expand_limits(y = c(0.5,1)) +
        theme(axis.text.y = element_text(size = 12)) +
        ##scale_y_continuous(label = function(x){paste0(x, ' (', as.integer(x/connections * 100), '%)')}) +
        ylab(paste0(description, ' (n = ', nodes, ')')) +
        xlab('False Positives, False Negatives (rates)') +
        geom_boxplot(outlier.size = 0) +
        geom_jitter(pch = 21, position = position_jitterdodge(jitter.height = 0.05, jitter.width = 2), alpha = .5) +
        ggtitle(title)
    
    return(p)
}

### end of file --

