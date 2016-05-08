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
    e_neg,
    nodes = NA,
    significance = 0.10,
    samples_significance = 0.001) {
    
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
    } else if (type == "random") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                random_dataset = sample.random.single.cells(i,e_pos[j],e_neg[j],nodes,significance,samples_significance)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = random_dataset$sampled_dataset
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
    wild_type,
    nodes = NA,
    significance = 0.10,
    samples_significance = 0.001) {
    
    
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
    } else if (type == "random") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                random_dataset = sample.random.multiple.biopses(i,e_pos[j],e_neg[j],nodes,significance,samples_significance,wild_type)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = random_dataset$sampled_dataset
                #results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
            }
        }
    }
    
    return(results)
}

# generate a random tree and the probabilities of the associated single cell samples
generate.random.single.cell.dataset <- function ( nodes, significance = 0.10 ) {
	
	# generate a random tree
	my_tree = generate.random.tree(nodes = nodes, significance = significance)
	
	# make an igraph object from the tree
	curr.igraph = graph.adjacency(my_tree$structure, mode = "directed")
	all_paths = get.all.shortest.paths(curr.igraph,my_tree$root,my_tree$leaves)$res
	my_paths = NULL
	for (i in 1:length(all_paths)) {
		my_paths[[i]] = as.vector(all_paths[[i]])
	}
	
	# build the dataset
	random_dataset = c(rep(0, nodes),(1-my_tree$root_prob))
	unique_random_dataset = c(rep(0, nodes),0.0)
	for (i in 1:length(my_paths)) {
		curr_path = my_paths[[i]]
		for (j in 1:length(curr_path)) {
			curr_valid_sample = NULL
			curr_sample_probability = 1
			for (k in 1:j) {
				curr_valid_sample = c(curr_valid_sample,curr_path[k])
				if(k==1) {
					curr_sample_probability = my_tree$root_prob
				}
				else {
					curr_sample_probability = curr_sample_probability * my_tree$probabilities[curr_path[k-1],curr_path[k]]
				}
			}
			new_valid_sample = rep(0,nodes)
			for(l in 1:length(curr_valid_sample)) {
				new_valid_sample[curr_valid_sample[l]] = 1
			}
			random_dataset = rbind(random_dataset,c(new_valid_sample,curr_sample_probability))
			unique_random_dataset = rbind(unique_random_dataset,c(new_valid_sample,0.0))
		}
	}
	unique_random_dataset = unique(unique_random_dataset)
	colnames(unique_random_dataset) = c(paste0("node_",1:(ncol(unique_random_dataset)-1)),"Probs")
	rownames(unique_random_dataset) = paste0("sample_",1:nrow(unique_random_dataset))
	
	# compute the probabilities for the unique random dataset
	for (i in 1:nrow(unique_random_dataset)) {
		curr_row = as.vector(unique_random_dataset[i,1:(ncol(unique_random_dataset)-1)])
		curr_prob = 0.0
		for (j in 1:nrow(random_dataset)) {
			curr_row_check = as.vector(random_dataset[j,1:(ncol(random_dataset)-1)])
			if(identical(curr_row,curr_row_check)) {
				curr_prob = curr_prob + random_dataset[j,ncol(random_dataset)]
			}
		}
		unique_random_dataset[i,ncol(unique_random_dataset)] = curr_prob
	}
	random_dataset = unique_random_dataset
	random_dataset[,"Probs"] = random_dataset[,"Probs"] / sum(random_dataset[,"Probs"])
	
	# set the names for the adjacency matrices
	colnames(my_tree$structure) = paste0("node_",1:ncol(my_tree$structure))
	rownames(my_tree$structure) = paste0("node_",1:ncol(my_tree$structure))
	colnames(my_tree$probabilities) = paste0("node_",1:ncol(my_tree$probabilities))
	rownames(my_tree$probabilities) = paste0("node_",1:ncol(my_tree$probabilities))
	
	# save the results
	my_random_dataset = list(structure = my_tree$structure, probabilities = my_tree$probabilities, 
	                            root_prob = my_tree$root_prob, dataset_samples = random_dataset)
	
	return(my_random_dataset)
	
}

# generate a random tree (both structure and conditional probabilities)
generate.random.tree <- function ( nodes, significance = 0.10 ) {
	
	# create the adjacency matrices to encode tree structure and probabilities
	adj.matrix.structure = array(0,c(nodes,nodes))
	adj.matrix.probabilities = array(0,c(nodes,nodes))
	
	# generate the tree structure randomly
	my_tree = generate.random.tree.structure(nodes)
	my_tree_structure = my_tree$structure
	my_tree_root = my_tree$root
	
	# create the actual structure of the tree
	root_prob = runif(1, min = significance, max = 1-significance)
	for (i in 2:length(my_tree_structure)) {
		curr_parents_nodes = my_tree_structure[[i-1]]
		curr_children_nodes = my_tree_structure[[i]]
		for (j in 1:length(curr_children_nodes)) {
			curr_node = curr_children_nodes[j]
			if(length(curr_parents_nodes)==1) {
				curr_parent = curr_parents_nodes[1]
			}
			else {
				curr_parent = sample(curr_parents_nodes,size=1)
			}
			curr_conditional_prob = runif(1, min = significance, max = 1-significance)
			# set the values in the adjacency matrices
			adj.matrix.structure[curr_parent,curr_node] = 1
			adj.matrix.probabilities[curr_parent,curr_node] = curr_conditional_prob
		}
	}
	
	# get the leaves in the tree
	leaves = NULL
	for (i in 1:nrow(adj.matrix.structure)) {
		nodes_children = which(adj.matrix.structure[i,]==1)
		if(length(nodes_children)==0) {
			leaves = c(leaves,i)
		}
	}
	
	# save the results
	my_tree = list(structure = adj.matrix.structure, probabilities = adj.matrix.probabilities, 
	                root = my_tree_root, leaves = leaves, root_prob = root_prob)
	
	return(my_tree)
	
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
	
	# save the results
	my_tree = list(structure = my_tree, root = my_tree[[1]])
	
	return(my_tree)
	
}
