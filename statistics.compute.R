##################################################################################
#                                                                                #
# MST                                                                            #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Giulio Caravagna, Luca De Sano, Daniele Ramazzotti         #
# email: tronco@disco.unimib.it                                                  #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################


# get the mean of the statistics
get.stats <- function(results) {

    experiments = colnames(results)
    sample_levels = rownames(results)
    noise_levels = names(results[[1,1]])
    my_algorithms = names(results[[1,1]][[1]][['reconstructions']])
    my_regularizators = c()

    for (algo in my_algorithms) {
        reg = names(results[[1,1]][[1]][['reconstructions']][[algo]])
        for (r in reg) {
            r = gsub('adj', 'res', r)

            my_regularizators = append(my_regularizators, r)
        }
    }
    my_regularizators = unique(my_regularizators)

    # structure to save the results
    my_results = NULL
    
    # get the statistics for the algorithms
    for (i in experiments) {
        for (a in my_algorithms) {
            for (b in my_regularizators) {


                if (!is.null(results[1,i][[1]][[as.character(1)]][["reconstructions"]][[a]][[b]])) {

                    cat('n ex: ', i, ' - algo: ', a, 'reg: ', b, '\n')
                    my.col.names = paste("Sample Level", sample_levels)
                    my.row.names = paste("Noise Level", noise_levels)

                    curr_result_accuracy = array(0,c(length(noise_levels),length(sample_levels)))
                    colnames(curr_result_accuracy) = my.col.names
                    rownames(curr_result_accuracy) = my.row.names
                    curr_result_sensitivity = array(0,c(length(noise_levels),length(sample_levels)))
                    colnames(curr_result_sensitivity) = my.col.names
                    rownames(curr_result_sensitivity) = my.row.names
                    curr_result_specificity = array(0,c(length(noise_levels),length(sample_levels)))
                    colnames(curr_result_specificity) = my.col.names
                    rownames(curr_result_specificity) = my.row.names
                    curr_result_hamming_distance = array(0,c(length(noise_levels),length(sample_levels)))
                    colnames(curr_result_hamming_distance) = my.col.names
                    rownames(curr_result_hamming_distance) = my.row.names
                    for (j in 1:length(sample_levels)) {
                        for (l in 1:length(noise_levels)) {
                            curr_computed = results[j,i][[1]][[as.character(l)]][["reconstructions"]][[a]][[b]]
                            curr_result_accuracy[l,j] = curr_computed[["accuracy"]]
                            curr_result_sensitivity[l,j] = curr_computed[["sensitivity"]]
                            curr_result_specificity[l,j] = curr_computed[["specificity"]]
                            curr_result_hamming_distance[l,j] = curr_computed[["hamming_distance"]]
                        }
                    }

                    if(i == 1) {
                        if (a == 'caprese') {
                            b = 'no.reg'
                        }

                        res = NULL
                        res[[1]] = list(curr_result_accuracy)
                        my_results[["accuracy"]][[a]][[gsub("\\.res","",b)]] = res
                        res = NULL
                        res[[1]] = list(curr_result_sensitivity)
                        my_results[["sensitivity"]][[a]][[gsub("\\.res","",b)]] = res
                        res = NULL
                        res[[1]] = list(curr_result_specificity)
                        my_results[["specificity"]][[a]][[gsub("\\.res","",b)]] = res
                        res = NULL
                        res[[1]] = list(curr_result_hamming_distance)
                        my_results[["hamming_distance"]][[a]][[gsub("\\.res","",b)]] = res
                        
                    } else {
                        if (a == 'caprese') {
                            b = 'no.reg'
                        }
                        
                        my_results[["accuracy"]][[a]][[gsub("\\.res","",b)]][[i]] = list(curr_result_accuracy)
                        my_results[["sensitivity"]][[a]][[gsub("\\.res","",b)]][[i]] = list(curr_result_sensitivity)
                        my_results[["specificity"]][[a]][[gsub("\\.res","",b)]][[i]] = list(curr_result_specificity)
                        my_results[["hamming_distance"]][[a]][[gsub("\\.res","",b)]][[i]] = list(curr_result_hamming_distance)
                        
                    }
                }

            }

#            if (a == "caprese") {
#                my.col.names = paste("Sample Level", sample_levels)
#                my.row.names = paste("Noise Level", noise_levels)
#                curr_result_accuracy = array(0,c(length(noise_levels),length(sample_levels)))
#                colnames(curr_result_accuracy) = my.col.names
#                rownames(curr_result_accuracy) = my.row.names
#                curr_result_sensitivity = array(0,c(length(noise_levels),length(sample_levels)))
#                colnames(curr_result_sensitivity) = my.col.names
#                rownames(curr_result_sensitivity) = my.row.names
#                curr_result_specificity = array(0,c(length(noise_levels),length(sample_levels)))
#                colnames(curr_result_specificity) = my.col.names
#                rownames(curr_result_specificity) = my.row.names
#                curr_result_hamming_distance = array(0,c(length(noise_levels),length(sample_levels)))
#                colnames(curr_result_hamming_distance) = my.col.names
#                rownames(curr_result_hamming_distance) = my.row.names
#
#
#
#
#
#                for (j in 1:length(sample_levels)) {
#                    for (l in 1:length(noise_levels)) {
#                        curr_computed = results[j,i][[1]][[as.character(l)]][["reconstructions"]][["caprese"]][["caprese.res"]]
#                        curr_result_accuracy[l,j] = curr_computed[["accuracy"]]
#                        curr_result_sensitivity[l,j] = curr_computed[["sensitivity"]]
#                        curr_result_specificity[l,j] = curr_computed[["specificity"]]
#                        curr_result_hamming_distance[l,j] = curr_computed[["hamming_distance"]]
#                    }
#                }
#                if (i == 1) {
#                    
#                    res = NULL
#                    res[[1]] = list(curr_result_accuracy)
#                    my_results[["accuracy"]][[a]][["no.reg"]] = res
#                    res = NULL
#                    res[[1]] = list(curr_result_sensitivity)
#                    my_results[["sensitivity"]][[a]][["no.reg"]] = res
#                    res = NULL
#                    res[[1]] = list(curr_result_specificity)
#                    my_results[["specificity"]][[a]][["no.reg"]] = res
#                    res = NULL
#                    res[[1]] = list(curr_result_hamming_distance)
#                    my_results[["hamming_distance"]][[a]][["no.reg"]] = res
#                    
#                } else {
#                        
#                    my_results[["accuracy"]][[a]][["no.reg"]][[i]] = list(curr_result_accuracy)
#                    my_results[["sensitivity"]][[a]][["no.reg"]][[i]] = list(curr_result_sensitivity)
#                    my_results[["specificity"]][[a]][["no.reg"]][[i]] = list(curr_result_specificity)
#                    my_results[["hamming_distance"]][[a]][["no.reg"]][[i]] = list(curr_result_hamming_distance)
#                    
#                }
#            }


        }
    }


    to.return = NULL
    for(stat.type in names(my_results)) {
        this.stat = my_results[[stat.type]]
        for (algorithm in names(this.stat)) {
            this.algorithm = this.stat[[algorithm]]
            for (reg in names(this.algorithm)) {
                this.reg = this.algorithm[[reg]]

                this.result = matrix(list(), nrow=length(noise_levels), ncol=length(sample_levels))
                this.stats = matrix(list(), nrow=length(noise_levels), ncol=length(sample_levels))

                for (experiment in 1:length(this.reg)) {
                    this.experiment = this.reg[[experiment]][[1]]
                    for (noise in 1:nrow(this.result)) {
                        for (sample in 1:ncol(this.result)) {
                            this.result[[noise, sample]] = c(this.result[[noise, sample]], this.experiment[[noise,sample]])
                        }
                    }
                }

                for (noise in 1:nrow(this.result)) {
                    for (sample in 1:ncol(this.result)) {
                        values = this.result[[noise, sample]]
                        mean = mean(values)
                        sd = sd(values)
                        median = median(values)
                        this.stats[[noise,sample]] = list(mean=mean, sd=sd, median=median, values=values)
                    }
                }

                to.return[[stat.type]][[algorithm]][[reg]] = this.stats
            }
        }
    }

    return(to.return)
}



