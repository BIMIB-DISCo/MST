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

# generate the structure of a random tree
generate.random.tree.structure <- function ( nodes ) {
	
	# if I have n nodes (>=2), I can have at minimum 2 levels (root plus level 2)
	# and at most an n-levels tree (when I have a path of n nodes)
	# select the number of levels randomly (uniform probability)
	if( nodes > 2 ) {
		all_nodes = 1:nodes
		levels = sample(2:nodes,size=1)
		my_tree = list()
		# put at least 1 element per level
		for (i in 1:levels) {
			if(length(all_nodes)>1) {
				curr_node = sample(all_nodes,size=1)
			}
			else {
				curr_node = all_nodes
			}
			my_tree[[i]] = curr_node
			all_nodes = all_nodes[-which(all_nodes==curr_node)]
		}
		# add the remaining edges
		while(length(all_nodes)>0) {
			if(length(all_nodes)>1) {
				curr_node = sample(all_nodes,size=1)
			}
			else {
				curr_node = all_nodes
			}
			if(levels==2) {
				curr_level = 2
			}
			else {
				curr_level = sample(2:levels,size=1)
			}
			my_tree[[curr_level]] = c(my_tree[[curr_level]],curr_node)
			all_nodes = all_nodes[-which(all_nodes==curr_node)]
		}
	}
	# in this case I have the minimum tree
	else if( nodes == 2 ) {
		my_tree = list()
		curr_nodes = sample(1:2,size=2)
		my_tree[[1]] = curr_nodes[1]
		my_tree[[2]] = curr_nodes[2]
	}
	# I need at least two nodes in my tree
	else {
		stop("I need at least 2 nodes in the tree!")
	}
	
	return(my_tree)
	
}
