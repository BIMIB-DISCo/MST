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


# generate a significant random tree with the Bayesian bulk model
sample.random.bayesian.bulk.significant.tree <- function ( nodes, min_significance = 0.6, max_significance = 0.90, sample_significance = 0.005 ) {
	
	# generate a random tree with minimum significance per samples
	random_tree = NULL
	while(is.null(random_tree) || 
	        min(random_tree$dataset_samples[,"Probs"])<=sample_significance || 
	        max(random_tree$dataset_samples[,"Probs"])>=(1-sample_significance)) {
	        	
		# keep generating random trees until I get a valid one
		random_tree = generate.random.bayesian.bulk.dataset(nodes,min_significance, max_significance)
		
	}
    
    return(random_tree)
	
}

# generate a random tree and the probabilities of the associated bulk cell samples
generate.random.bayesian.bulk.dataset <- function ( nodes, min_significance = 0.6, max_significance = 0.9 ) {
    
    # generate a random tree
    my_tree = generate.random.tree(nodes = nodes, min_significance = min_significance, max_significance = max_significance)
    
    # get adj.matrix
    adj.matrix = my_tree$structure
    root = which(colSums(adj.matrix) == 0)
    
    # get probabilities
    probabilities = my_tree$probabilities
    root.prob = my_tree$root_prob
    
    # get the parent of each node
    parents.pos = rep(-1,nodes)
    for(n in 1:nodes) {
        if(n!=root) {
            parents.pos[n] = which(probabilities[,n]>0)
        }
    }

    # Build first row of sample  datasera (1 - p(root))
    random_dataset = c(rep(0, nodes),(1-my_tree$root_prob))
    
    # consider all the possible binary strings of #nodes bits
    my.data = lapply(1:(2^nodes-1),FUN=function(i) {
        # get a valid sample
        is.valid = TRUE
        curr_sample = c(rep(0,nodes),0.0)
        curr_binary = binary(i)
        curr_sample[(length(curr_sample)-length(curr_binary)):(length(curr_sample)-1)] = curr_binary
        is.valid = check.valid(curr_sample,parents.pos)
        # if the sample is valid, compute its probability and add it to the dataset
        if(is.valid==TRUE) {
            curr_prob = 1
            for(j in which(curr_sample==1)) {
                if(j==root) {
                    curr_prob = curr_prob * root.prob
                }
                else {
                    curr_prob = curr_prob * probabilities[parents.pos[j],j]
                }
            }
            curr_sample[length(curr_sample)] = curr_prob
        }
        else {
            curr_sample = NULL
        }
        return(curr_sample)
    })
    my.data = Reduce("rbind",my.data)
    random_dataset = rbind(random_dataset,my.data)
    colnames(random_dataset) = c(paste0("node_",1:(ncol(random_dataset)-1)),"Probs")
    rownames(random_dataset) = paste0("sample_",1:nrow(random_dataset))

    # normalize the probabilities
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

# convert from decimal to binary
binary <- function(n) {
    binary = NULL
    if(n > 1) {
        binary = c(binary,binary(as.integer(n/2)))
    }
    binary = c(binary,n %% 2)
    return(binary)
}

# check if a sample is valid according to the generative model
check.valid = function( sample, parents.pos ) {
    is.valid = TRUE
    is.valid = all(sapply(which(sample==1),FUN=function( x ) {
        curr.valid = TRUE
        if(parents.pos[x]!=-1) {
            if(sample[parents.pos[x]]==0) {
                curr.valid = FALSE
            }
        }
        return(curr.valid)
    }))
    return(is.valid)
}
