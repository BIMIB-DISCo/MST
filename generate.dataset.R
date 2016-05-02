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


# simulate a dataset from single cells at given sample sizes and noise levels
generate.dataset.single.cells <- function (type,
    true_tree,
    samples_num,
    nodes_probabilities = NA,
    e_pos,
    e_neg) {
    
    # structure to save the results
    results = NULL
    
    # run the experiments
    if (type == "low") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.single.cells.polyclonal.low(i,nodes_probabilities,e_pos[j],e_neg[j])
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                #results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
            }
        }
    } else if (type == "medium") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.single.cells.polyclonal.medium(i,nodes_probabilities,e_pos[j],e_neg[j])
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                #results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
            }
        }
    } else if (type == "high") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.single.cells.polyclonal.high(i,nodes_probabilities,e_pos[j],e_neg[j])
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                #results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
            }
        }
    }
    
    return(results)
    
}


# simulate a dataset from multiple biopses at given sample sizes and noise levels
generate.dataset.multiple.biopses <- function(type,
    true_tree,
    samples_num,
    clones_probabilities,
    nodes_probabilities = NA,
    e_pos,
    e_neg,
    wild_type) {
    
    
    # structure to save the results
    results = NULL
    
    # run the experiments
    if (type == "low") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.multiple.biopses.polyclonal.low(i,
                    clones_probabilities,
                    nodes_probabilities,
                    e_pos[j],
                    e_neg[j],
                    wild_type)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                #results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
            }
        }
    } else if (type == "medium") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.multiple.biopses.polyclonal.medium(i,
                    clones_probabilities,
                    nodes_probabilities,
                    e_pos[j],
                    e_neg[j],
                    wild_type)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                #results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
            }
        }
    } else if (type == "high") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.multiple.biopses.polyclonal.high(i,
                    clones_probabilities,
                    nodes_probabilities,
                    e_pos[j],
                    e_neg[j],
                    wild_type)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                #results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
            }
        }
    }
    return(results)
}
