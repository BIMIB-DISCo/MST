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
    true_tree = NULL,
    samples_num,
    nodes_probabilities = NA,
    e_pos,
    e_neg,
    nodes = NA,
    min_significance = 0.6,
    max_significance = 0.9,
    samples_significance = 0.001) {
        
    
    library(igraph)
    
    # structure to save the results
    results = NULL
    
    # run the experiments
    if (type == "low") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.single.cells.polyclonal.low(i,nodes_probabilities,e_pos[j],e_neg[j])
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = true_tree
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]

            }
        }
    } else if (type == "medium") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.single.cells.polyclonal.medium(i,nodes_probabilities,e_pos[j],e_neg[j])
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = true_tree
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
            }
        }
    } else if (type == "high") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                curr_dataset = sample.single.cells.polyclonal.high(i,nodes_probabilities,e_pos[j],e_neg[j])
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = true_tree
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
            }
        }
    } else if (type == "random") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                random_dataset = sample.random.single.cells(i,e_pos[j],e_neg[j],nodes,min_significance,max_significance,samples_significance)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = random_dataset$sampled_dataset
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = random_dataset$random_tree$structure
                results[[as.character(i)]][[as.character(j)]][["probabilities"]] = random_dataset$random_tree$probabilities
                results[[as.character(i)]][[as.character(j)]][["root_prob"]] = random_dataset$random_tree$root_prob
                results[[as.character(i)]][[as.character(j)]][["dataset_samples"]] = random_dataset$random_tree$dataset_samples
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
            }
        }
    } else if (type == "random_forest") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                random_dataset = sample.random.single.cells.forest(i, e_pos[j],e_neg[j],nodes,min_significance,max_significance,samples_significance)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = random_dataset$sampled_dataset
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = random_dataset$random_tree$structure
                results[[as.character(i)]][[as.character(j)]][["dataset_samples"]] = random_dataset$random_tree$dataset_samples
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
            }
        }
    }
    
    return(results)
    
}


# simulate a dataset from multiple biopses at given sample sizes and noise levels
generate.dataset.multiple.biopses <- function(type,
    true_tree = NULL,
    samples_num,
    clones_probabilities,
    nodes_probabilities = NA,
    e_pos,
    e_neg,
    wild_type,
    nodes = NA,
    min_significance = 0.60,
    max_significance = 0.90,
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
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = true_tree
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
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
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = true_tree
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
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
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = true_tree
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
            }
        }
    } else if (type == "random") {
        for (i in samples_num) {
            for (j in 1:length(e_pos)) {
                random_dataset = sample.random.multiple.biopses(i,e_pos[j],e_neg[j],nodes,min_significance,max_significance,samples_significance,wild_type)
                results[[as.character(i)]][[as.character(j)]][["dataset"]] = random_dataset$sampled_dataset
                results[[as.character(i)]][[as.character(j)]][["true_tree"]] = random_dataset$random_tree$structure
                results[[as.character(i)]][[as.character(j)]][["dataset_samples"]] = random_dataset$random_tree$dataset_samples
                results[[as.character(i)]][[as.character(j)]][["epos"]] = e_pos[j]
                results[[as.character(i)]][[as.character(j)]][["eneg"]] = e_neg[j]
            }
        }
    }
    
    return(results)
}

# generate a random tree and the probabilities of the associated single cell samples
generate.random.single.cell.dataset <- function ( nodes, min_significance = 0.6, max_significance = 0.90 ) {
    
    # generate a random tree
    my_tree = generate.random.tree(nodes = nodes, min_significance = min_significance, max_significance = max_significance)
    
    # get adj.matrix
    adj.matrix = my_tree$structure
    root = which(colSums(adj.matrix) == 0)
    leaves = which(rowSums(adj.matrix) == 0)
    
    # make an igraph object from the tree
    graph = graph.adjacency(my_tree$structure, mode = "directed")   
    paths = lapply(get.all.shortest.paths(graph,root,leaves)$res, as.vector)
    
    # get probabilities
    probabilities = my_tree$probabilities
    root.prob = my_tree$root_prob

    # evaluate the probability of not having children
    not.having.children <- function(node) {
        children = which(adj.matrix[node,] == 1)
        children.prob  = 1
        for (child in children) {
            children.prob = children.prob * (1 - probabilities[node, child])
        }
        return(children.prob)
    }

    # evaluate the probability of not having brothers
    not.having.brothers <- function(node, parent) {
        brothers = which(adj.matrix[parent,] == 1)
        brothers = brothers[brothers != node]
        brothers.prob  = 1
        for (brother in brothers) {
            brothers.prob = brothers.prob * (1 - probabilities[parent, brother])
        }
        return(brothers.prob)
    }

    # if I move from the root, the probability of the sample is 
    # P(Root) * P(Sample) * (1 - P(alternative_paths)), where
    # alternative_paths are alternatives to the chosen clone

    # Build first row of sample  datasera (1 - p(root))
    random_dataset = c(rep(0, nodes),(1-my_tree$root_prob))
    
    # consider all the possible clones described by the tree
    for (path in paths) {
        # go through all the possible samples in this clonal path
        for (j in 1:length(path)) {
            sample = NULL
            sample_probability = 1
            # visit any subset of the current clonal path
            for (k in 1:j) {
                sample = c(sample, path[k])
                # Now check every possible configuration

                # If this node is the root ... 
                if(path[k] == root) {
                    # If there is only this node...
                    if (k == j) {
                        # compute the probabilities of not having
                        # children event
                        children.prob = not.having.children(path[k])
                        sample_probability = root.prob * children.prob
                    } else {
                        # ...else this is the first step of a progression
                        sample_probability = root.prob
                    }

                # ...else if this node is a leaf...
                } else if (path[k] %in% leaves) {
                    # compute the probabilities of not having
                    # brother event
                    brothers.prob = not.having.brothers(path[k], path[k-1])
                    sample_probability = sample_probability * probabilities[path[k-1], path[k]] * brothers.prob
                
                # ...else this node is in the middle of a progression
                } else {
                    # If this is the last node of a subprogression
                    if (k == j) {
                        #print('this path terminate here')
                        brothers.prob = not.having.brothers(path[k], path[k-1])
                        children.prob = not.having.children(path[k])
                        sample_probability = sample_probability * probabilities[path[k-1], path[k]] * brothers.prob * children.prob
                    } else {
                        #print('we are in the middle of a path')
                        brothers.prob = not.having.brothers(path[k], path[k-1])
                        sample_probability = sample_probability * probabilities[path[k-1], path[k]] * brothers.prob
                    }
                }
            }
                            
            valid_sample = rep(0,nodes)
            for(l in 1:length(sample)) {
                valid_sample[sample[l]] = 1
            }
            random_dataset = rbind(random_dataset, c(valid_sample,sample_probability))
        }
    }
    random_dataset = unique(random_dataset)
    colnames(random_dataset) = c(paste0("node_",1:(ncol(random_dataset)-1)),"Probs")
    rownames(random_dataset) = paste0("sample_",1:nrow(random_dataset))

    # normalize a la eros
    random_dataset[2:nrow(random_dataset),"Probs"] = random_dataset[2:nrow(random_dataset),"Probs"] / sum(random_dataset[2:nrow(random_dataset),"Probs"]) * my_tree$root_prob

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
generate.random.tree <- function ( nodes, min_significance = 0.60, max_significance = 0.90 ) {
    
    # create the adjacency matrices to encode tree structure and probabilities
    adj.matrix.structure = array(0,c(nodes,nodes))
    adj.matrix.probabilities = array(0,c(nodes,nodes))
    
    # generate the tree structure randomly
    my_tree = generate.random.tree.structure(nodes)
    my_tree_structure = my_tree$structure
    my_tree_root = my_tree$root
    
    # create the actual structure of the tree
    root_prob = runif(1, min = min_significance, max = max_significance)
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
            curr_conditional_prob = runif(1, min = min_significance, max = max_significance)
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
