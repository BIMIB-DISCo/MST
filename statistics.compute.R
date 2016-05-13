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
get.stats <- function(results,
	my_algorithms,
	my_regularizators,
	sample_levels,
	noise_levels,
	number_experiments) {
	
	# structure to save the results
	my_results = NULL
	
	# get the statistics for the algorithms
	for (i in 1:number_experiments) {
		for (a in my_algorithms) {
			for (b in my_regularizators) {
				if (!is.null(results[1,i][[1]][[as.character(1)]][["reconstructions"]][[a]][[b]])) {
					my.col.names = paste("Sample Level",1:sample_levels)
					my.row.names = paste("Noise Level",1: noise_levels)
					curr_result_accuracy = array(0,c(noise_levels,sample_levels))
					colnames(curr_result_accuracy) = my.col.names
					rownames(curr_result_accuracy) = my.row.names
					curr_result_sensitivity = array(0,c(noise_levels,sample_levels))
					colnames(curr_result_sensitivity) = my.col.names
					rownames(curr_result_sensitivity) = my.row.names
					curr_result_specificity = array(0,c(noise_levels,sample_levels))
					colnames(curr_result_specificity) = my.col.names
					rownames(curr_result_specificity) = my.row.names
					for (j in 1:sample_levels) {
						for (l in 1:noise_levels) {
							curr_computed = results[j,i][[1]][[as.character(l)]][["reconstructions"]][[a]][[b]]
							curr_result_accuracy[l,j] = curr_computed[["accuracy"]]
							curr_result_sensitivity[l,j] = curr_computed[["sensitivity"]]
							curr_result_specificity[l,j] = curr_computed[["specificity"]]
							curr_result_hamming_distance[l,j] = curr_computed[["hamming_distance"]]
						}
					}
					if(i == 1) {
						
						res = NULL
						res[[1]] = list(curr_result_accuracy)
						my_results[["accuracy"]][[a]][[gsub(".res","",b)]] = res
						res = NULL
						res[[1]] = list(curr_result_sensitivity)
						my_results[["sensitivity"]][[a]][[gsub(".res","",b)]] = res
						res = NULL
						res[[1]] = list(curr_result_specificity)
						my_results[["specificity"]][[a]][[gsub(".res","",b)]] = res
						res = NULL
						res[[1]] = list(curr_result_hamming_distance)
						my_results[["hamming_distance"]][[a]][[gsub(".res","",b)]] = res
						
					} else {
						
						my_results[["accuracy"]][[a]][[gsub(".res","",b)]][[i]] = list(curr_result_accuracy)
						my_results[["sensitivity"]][[a]][[gsub(".res","",b)]][[i]] = list(curr_result_sensitivity)
						my_results[["specificity"]][[a]][[gsub(".res","",b)]][[i]] = list(curr_result_specificity)
						my_results[["hamming_distance"]][[a]][[gsub(".res","",b)]][[i]] = list(curr_result_hamming_distance)
						
					}
				}
			}
			if (a == "caprese") {
				my.col.names = paste("Sample Level",1:sample_levels)
				my.row.names = paste("Noise Level",1: noise_levels)
				curr_result_accuracy = array(0,c(noise_levels,sample_levels))
				colnames(curr_result_accuracy) = my.col.names
				rownames(curr_result_accuracy) = my.row.names
				curr_result_sensitivity = array(0,c(noise_levels,sample_levels))
				colnames(curr_result_sensitivity) = my.col.names
				rownames(curr_result_sensitivity) = my.row.names
				curr_result_specificity = array(0,c(noise_levels,sample_levels))
				colnames(curr_result_specificity) = my.col.names
				rownames(curr_result_specificity) = my.row.names
				for (j in 1:sample_levels) {
					for (l in 1:noise_levels) {
						curr_computed = results[j,i][[1]][[as.character(l)]][["reconstructions"]][["caprese"]][["caprese.res"]]
						curr_result_accuracy[l,j] = curr_computed[["accuracy"]]
						curr_result_sensitivity[l,j] = curr_computed[["sensitivity"]]
						curr_result_specificity[l,j] = curr_computed[["specificity"]]
						curr_result_hamming_distance[l,j] = curr_computed[["hamming_distance"]]
					}
				}
				if (i == 1) {
					
					res = NULL
					res[[1]] = list(curr_result_accuracy)
					my_results[["accuracy"]][[a]][["no.reg"]] = res
					res = NULL
					res[[1]] = list(curr_result_sensitivity)
					my_results[["sensitivity"]][[a]][["no.reg"]] = res
					res = NULL
					res[[1]] = list(curr_result_specificity)
					my_results[["specificity"]][[a]][["no.reg"]] = res
					res = NULL
					res[[1]] = list(curr_result_hamming_distance)
					my_results[["hamming_distance"]][[a]][["no.reg"]] = res
					
				} else {
						
					my_results[["accuracy"]][[a]][["no.reg"]][[i]] = list(curr_result_accuracy)
					my_results[["sensitivity"]][[a]][["no.reg"]][[i]] = list(curr_result_sensitivity)
					my_results[["specificity"]][[a]][["no.reg"]][[i]] = list(curr_result_specificity)
					my_results[["hamming_distance"]][[a]][["no.reg"]][[i]] = list(curr_result_hamming_distance)
					
				}
			}
		}
	}
	
	# mean the results
	for (a in names(my_results)) {
		for (b in names(my_results[[a]])) {
			for (c in names(my_results[[a]][[b]])) {
				curr.res = my_results[[a]][[b]][[c]]
				curr_my_res = array(list(),c(nrow(curr.res[[1]][[1]]),ncol(curr.res[[1]][[1]])))
				# consider all the combinations sample size / noise level
				for (i in 1:nrow(curr.res[[1]][[1]])) {
					for(j in 1:ncol(curr.res[[1]][[1]])) {
						
						# consider all the different experiments
						curr_vals = NULL
						for (all_curr_exp in 1:length(curr.res)) {
							curr_vals = c(curr_vals,curr.res[[all_curr_exp]][[1]][i,j])
						}
						
						# compute the stats
						curr.res_mean = mean(curr_vals)
						curr.res_sd = sd(curr_vals)
						curr_my_res[[i,j]] = list(mean=curr.res_mean,sd=curr.res_sd)
						
					}
				}
				my_results[[a]][[b]][[c]] = curr_my_res
			}
		}
	}
	return(my_results)
}
