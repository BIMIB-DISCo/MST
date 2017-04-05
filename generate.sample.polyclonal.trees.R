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

# sample single cells from the polyclonal-low tree model
sample.single.cells.polyclonal.low <- function (samples_num,
    nodes_probabilities = NA,
    e_pos,
    e_neg) {
    
    # create the structure with the valid samples
    dataset = matrix(0, 7, 6)
    dataset[2,1] = 1
    dataset[3,1] = 1
    dataset[3,2] = 1
    dataset[4,1] = 1
    dataset[4,3] = 1
    dataset[5,1] = 1
    dataset[5,4] = 1
    dataset[6,1] = 1
    dataset[6,2] = 1
    dataset[6,5] = 1
    dataset[7,1] = 1
    dataset[7,3] = 1
    dataset[7,6] = 1
    
    # randomply sample the nodes probabilities if not provided
    if (length(nodes_probabilities) == 0 || is.na(nodes_probabilities)) {
        nodes_probabilities = runif(6, 0.4, 0.6)
    }
    
    # set the probabilities of the nodes
    p1 = nodes_probabilities[1]
    p2 = nodes_probabilities[2]
    p3 = nodes_probabilities[3]
    p4 = nodes_probabilities[4]
    p5 = nodes_probabilities[5]
    p6 = nodes_probabilities[6]
    
    # set the probabilities of the samples
    samples_probabilities = rep(0,7)
    samples_probabilities[1] = 1 - p1
    samples_probabilities[2] = p1 * (1 - p2) * (1 - p3) * (1 - p4)
    samples_probabilities[3] = p1 * p2 * (1 - p3) * (1 - p4) * (1 - p5)
    samples_probabilities[4] = p1 * (1 - p2) * p3 * (1 - p4) * (1 - p6)
    samples_probabilities[5] = p1 * (1 - p2) * (1 - p3) * p4
    samples_probabilities[6] = p1 * p2 * (1 - p3) * (1 - p4) * p5
    samples_probabilities[7] = p1 * (1 - p2) * p3 * (1 - p4) * p6
    
    # normalize the samples probabilities
    samples_probabilities = samples_probabilities / sum(samples_probabilities)
    
    # sample the dataset
    sampled_dataset = dataset[sample(1:nrow(dataset),
        size = samples_num,
        replace = TRUE,
        prob = samples_probabilities),]
    
    # apply the noise
    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,],e_pos,e_neg)
    }
    
    return(sampled_dataset)
    
}

# sample single cells from the polyclonal-medium tree model
sample.single.cells.polyclonal.medium <- function (samples_num,
    nodes_probabilities = NA,
    e_pos,
    e_neg) {
    
    # create the structure with the valid samples
    dataset = matrix(0, 12, 11)
    dataset[2,1] = 1
    dataset[3,1] = 1
    dataset[3,2] = 1
    dataset[4,1] = 1
    dataset[4,3] = 1
    dataset[5,1] = 1
    dataset[5,4] = 1
    dataset[6,1] = 1
    dataset[6,2] = 1
    dataset[6,5] = 1
    dataset[7,1] = 1
    dataset[7,4] = 1
    dataset[7,6] = 1
    dataset[8,1] = 1
    dataset[8,2] = 1
    dataset[8,5] = 1
    dataset[8,7] = 1
    dataset[9,1] = 1
    dataset[9,2] = 1
    dataset[9,5] = 1
    dataset[9,8] = 1
    dataset[10,1] = 1
    dataset[10,4] = 1
    dataset[10,6] = 1
    dataset[10,9] = 1
    dataset[11,1] = 1
    dataset[11,4] = 1
    dataset[11,6] = 1
    dataset[11,10] = 1
    dataset[12,1] = 1
    dataset[12,4] = 1
    dataset[12,6] = 1
    dataset[12,11] = 1
    
    # randomply sample the nodes probabilities if not provided
    if (length(nodes_probabilities) == 0 || is.na(nodes_probabilities)) {
        nodes_probabilities = runif(11, 0.4, 0.6)
    }
    
    # set the probabilities of the nodes
    p1 = nodes_probabilities[1]
    p2 = nodes_probabilities[2]
    p3 = nodes_probabilities[3]
    p4 = nodes_probabilities[4]
    p5 = nodes_probabilities[5]
    p6 = nodes_probabilities[6]
    p7 = nodes_probabilities[7]
    p8 = nodes_probabilities[8]
    p9 = nodes_probabilities[9]
    p10 = nodes_probabilities[10]
    p11 = nodes_probabilities[11]
    
    # set the probabilities of the samples
    samples_probabilities = rep(0,12)
    samples_probabilities[1] = 1 - p1
    samples_probabilities[2] = p1 * (1 - p2) * (1 - p3) * (1 - p4)
    samples_probabilities[3] = p1 * p2 * (1 - p3) * (1 - p4) * (1 - p5)
    samples_probabilities[4] = p1 * (1 - p2) * p3 * (1 - p4)
    samples_probabilities[5] = p1 * (1 - p2) * (1 - p3) * p4 * (1 - p6)
    samples_probabilities[6] = p1 * p2 * (1 - p3) * (1 - p4) * p5 * (1 - p7) * (1 - p8)
    samples_probabilities[7] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * (1 - p10) * (1 - p11)
    samples_probabilities[8] = p1 * p2 * (1 - p3) * (1 - p4) * p5 * p7 * (1 - p8)
    samples_probabilities[9] = p1 * p2 * (1 - p3) * (1 - p4) * p5 * (1 - p7) * p8
    samples_probabilities[10] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * p9 * (1 - p10) * (1 - p11)
    samples_probabilities[11] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * p10 * (1 - p11)
    samples_probabilities[12] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * (1 - p10) * p11
    
    # normalize the samples probabilities
    samples_probabilities = samples_probabilities / sum(samples_probabilities)
    
    # sample the dataset
    sampled_dataset = dataset[sample(1:nrow(dataset),
        size = samples_num,
        replace = TRUE,
        prob = samples_probabilities),]
    
    # apply the noise
    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,],e_pos,e_neg)
    }
    
    return(sampled_dataset)
    
}

# sample single cells from the polyclonal-high tree model
sample.single.cells.polyclonal.high <- function (samples_num,
    nodes_probabilities = NA,
    e_pos,
    e_neg) {
    
    # create the structure with the valid samples
    dataset = matrix(0, 18, 17)
    dataset[2,1] = 1
    dataset[3,1] = 1
    dataset[3,2] = 1
    dataset[4,1] = 1
    dataset[4,3] = 1
    dataset[5,1] = 1
    dataset[5,4] = 1
    dataset[6,1] = 1
    dataset[6,2] = 1
    dataset[6,5] = 1
    dataset[7,1] = 1
    dataset[7,4] = 1
    dataset[7,6] = 1
    dataset[8,1] = 1
    dataset[8,2] = 1
    dataset[8,5] = 1
    dataset[8,7] = 1
    dataset[9,1] = 1
    dataset[9,2] = 1
    dataset[9,5] = 1
    dataset[9,8] = 1
    dataset[10,1] = 1
    dataset[10,4] = 1
    dataset[10,6] = 1
    dataset[10,9] = 1
    dataset[11,1] = 1
    dataset[11,4] = 1
    dataset[11,6] = 1
    dataset[11,10] = 1
    dataset[12,1] = 1
    dataset[12,4] = 1
    dataset[12,6] = 1
    dataset[12,11] = 1
    dataset[13,1] = 1
    dataset[13,2] = 1
    dataset[13,5] = 1
    dataset[13,7] = 1
    dataset[13,12] = 1
    dataset[14,1] = 1
    dataset[14,4] = 1
    dataset[14,6] = 1
    dataset[14,9] = 1
    dataset[14,13] = 1
    dataset[15,1] = 1
    dataset[15,4] = 1
    dataset[15,6] = 1
    dataset[15,9] = 1
    dataset[15,14] = 1
    dataset[16,1] = 1
    dataset[16,4] = 1
    dataset[16,6] = 1
    dataset[16,9] = 1
    dataset[16,15] = 1
    dataset[17,1] = 1
    dataset[17,4] = 1
    dataset[17,6] = 1
    dataset[17,11] = 1
    dataset[17,16] = 1
    dataset[18,1] = 1
    dataset[18,4] = 1
    dataset[18,6] = 1
    dataset[18,11] = 1
    dataset[18,17] = 1
    
    # randomply sample the nodes probabilities if not provided
    if (length(nodes_probabilities) == 0 || is.na(nodes_probabilities)) {
        nodes_probabilities = runif(17, 0.4, 0.6)
    }
    
    # set the probabilities of the nodes
    p1 = nodes_probabilities[1]
    p2 = nodes_probabilities[2]
    p3 = nodes_probabilities[3]
    p4 = nodes_probabilities[4]
    p5 = nodes_probabilities[5]
    p6 = nodes_probabilities[6]
    p7 = nodes_probabilities[7]
    p8 = nodes_probabilities[8]
    p9 = nodes_probabilities[9]
    p10 = nodes_probabilities[10]
    p11 = nodes_probabilities[11]
    p12 = nodes_probabilities[12]
    p13 = nodes_probabilities[13]
    p14 = nodes_probabilities[14]
    p15 = nodes_probabilities[15]
    p16 = nodes_probabilities[16]
    p17 = nodes_probabilities[17]
    
    # set the probabilities of the samples
    samples_probabilities = rep(0, 18)
    samples_probabilities[1] = 1 - p1
    samples_probabilities[2] = p1 * (1 - p2) * (1 - p3) * (1 - p4)
    samples_probabilities[3] = p1 * p2 * (1 - p3) * (1 - p4) * (1 - p5)
    samples_probabilities[4] = p1 * (1 - p2) * p3 * (1 - p4)
    samples_probabilities[5] = p1 * (1 - p2) * (1 - p3) * p4 * (1 - p6)
    samples_probabilities[6] = p1 * p2 * (1 - p3) * (1 - p4) * p5 * (1 - p7) * (1 - p8)
    samples_probabilities[7] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * (1 - p10) * (1 - p11)
    samples_probabilities[8] = p1 * p2 * (1 - p3) * (1 - p4) * p5 * p7 * (1 - p8) * (1 - p12)
    samples_probabilities[9] = p1 * p2 * (1 - p3) * (1 - p4) * p5 * (1 - p7) * p8
    samples_probabilities[10] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * p9 * (1 - p10) * (1 - p11) * (1 - p13) * (1 - p14) * (1 - p15)
    samples_probabilities[11] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * p10 * (1 - p11)
    samples_probabilities[12] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * (1 - p10) * p11 * (1 - p16) * (1 - p17)
    samples_probabilities[13] = p1 * p2 * (1 - p3) * (1 - p4) * p5 * p7 * (1 - p8) * p12
    samples_probabilities[14] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * p9 * (1 - p10) * (1 - p11) * p13 * (1 - p14) * (1 - p15)
    samples_probabilities[15] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * p9 * (1 - p10) * (1 - p11) * (1 - p13) * p14 * (1 - p15)
    samples_probabilities[16] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * p9 * (1 - p10) * (1 - p11) * (1 - p13) * (1 - p14) * p15
    samples_probabilities[17] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * (1 - p10) * p11 * p16 * (1 - p17)
    samples_probabilities[18] = p1 * (1 - p2) * (1 - p3) * p4 * p6 * (1 - p9) * (1 - p10) * p11 * (1 - p16) * p17
    
    # normalize the samples probabilities
    samples_probabilities = samples_probabilities / sum(samples_probabilities)
    
    # sample the dataset
    sampled_dataset = dataset[sample(1:nrow(dataset),
        size = samples_num,
        replace = TRUE,
        prob = samples_probabilities),]
    
    # apply the noise
    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }
    
    return(sampled_dataset)
    
}

# generate a dataset of random single cell trees
sample.random.single.cells <- function ( samples_num, e_pos, e_neg, nodes, min_significance = 0.60, max_significance = 0.9, sample_significance = 0.005, true_tree = NULL ) {
	# generate a random tree with minimum significance per samples
	random_tree = NULL
	while(is.null(random_tree) || 
	        min(random_tree$dataset_samples[,"Probs"])<=sample_significance || 
	        max(random_tree$dataset_samples[,"Probs"])>=(1-sample_significance)) {
	        	

		# keep generating random trees until I get a valid one
		random_tree = generate.random.single.cell.dataset(nodes,min_significance,max_significance,true_tree)
		
	}
	
	dataset = random_tree$dataset_samples[,1:nodes]
	samples_probabilities = random_tree$dataset_samples[,nodes+1]
	
	# sample the dataset
    sampled_dataset = dataset[sample(1:nrow(dataset),
        size = samples_num,
        replace = TRUE,
        prob = samples_probabilities),]
    
    # apply the noise
    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }
    rownames(sampled_dataset) = paste0("sample_",1:samples_num)
    
    # save the results
    res = list( random_tree = random_tree, sampled_dataset = sampled_dataset )
    
    return(res)
	
}

# generate a significant random tree
sample.random.single.cells.significant.tree <- function ( nodes, min_significance = 0.6, max_significance = 0.90, sample_significance = 0.005, true_tree = NULL ) {
	
	# generate a random tree with minimum significance per samples
	random_tree = NULL
	while(is.null(random_tree) || 
	        min(random_tree$dataset_samples[,"Probs"])<=sample_significance || 
	        max(random_tree$dataset_samples[,"Probs"])>=(1-sample_significance)) {
	        	
		# keep generating random trees until I get a valid one
		random_tree = generate.random.single.cell.dataset(nodes,min_significance, max_significance, true_tree)
		
	}
    
    return(random_tree)
	
}

# generate a dataset of random single cell given the generating tree
sample.random.single.cells.given.tree <- function ( samples_num, nodes, e_pos, e_neg, random_tree ) {
	
	dataset = random_tree$dataset_samples[,1:nodes]
	samples_probabilities = random_tree$dataset_samples[,nodes+1]
	
	# sample the dataset
    sampled_dataset = dataset[sample(1:nrow(dataset),
        size = samples_num,
        replace = TRUE,
        prob = samples_probabilities),]
    
    # apply the noise
    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }
    rownames(sampled_dataset) = paste0("sample_",1:samples_num)
    
    sampled_dataset = colSums(sampled_dataset) / samples_num
    
    # if a mutation is observed at least in 1 sample, it is observed (OR of the samples)
    sampled_dataset[which(sampled_dataset > 0)] = 1
    
    return(sampled_dataset)
	
}

# sample one observation of multiple biospes from a given single cell model
sample.multiple.biospes.from.single.cells <- function (clones_num,
    single_cell_type,
    nodes_probabilities = NA,
    e_pos,
    e_neg) {
    
    # generate the single cells data to be combined
    if (single_cell_type == "low") {
        sampled_dataset = sample.single.cells.polyclonal.low(clones_num, nodes_probabilities, e_pos, e_neg)
    } else if (single_cell_type == "medium") {
        sampled_dataset = sample.single.cells.polyclonal.medium(clones_num, nodes_probabilities, e_pos, e_neg)
    } else if (single_cell_type == "high") {
        sampled_dataset = sample.single.cells.polyclonal.high(clones_num, nodes_probabilities, e_pos, e_neg)
    }
    
    sampled_dataset = colSums(sampled_dataset) / clones_num
    
    # if a mutation is observed at least in 1 sample, it is observed (OR of the samples)
    sampled_dataset[which(sampled_dataset > 0)] = 1
    
    # return(colSums(sampled_dataset))
    return(sampled_dataset)
}

# sample multiple biopses from the polyclonal-low tree model
sample.multiple.biopses.polyclonal.low <- function (samples_num,
    clones_probabilities,
    nodes_probabilities = NA,
    e_pos,
    e_neg,
    wild_type) {
    
    wild_type_samples = round(samples_num * wild_type)
    if (wild_type > 0 && wild_type_samples == 0) {
        wild_type_samples = 1
    }
    res.matrix = matrix(NA, ncol=6, nrow=(samples_num - wild_type_samples))
    sampled_dataset = sapply(res.matrix[,1], function(x){
        sample.multiple.biospes.from.single.cells(clones_probabilities,
            "low",
            nodes_probabilities,
            0,
            0)
    })
    sampled_dataset = t(sampled_dataset)
    
    if (wild_type_samples > 0) {
        sampled_dataset = rbind(matrix(0, wild_type_samples,6), sampled_dataset)
    }

    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }

    return(sampled_dataset)
}

# sample multiple biopses from the polyclonal-medium tree model
sample.multiple.biopses.polyclonal.medium <- function(samples_num,
    clones_probabilities,
    nodes_probabilities = NA,
    e_pos,
    e_neg,
    wild_type) {
    
    wild_type_samples = round(samples_num * wild_type)
    if (wild_type > 0 && wild_type_samples == 0) {
        wild_type_samples = 1
    }
    
    res.matrix = matrix(NA, ncol=11, nrow=(samples_num-wild_type_samples))
    sampled_dataset = sapply(res.matrix[,1], function(x){
        sample.multiple.biospes.from.single.cells(clones_probabilities,
            "medium",
            nodes_probabilities,
            0,
            0)
    })
    sampled_dataset = t(sampled_dataset)
    
    if (wild_type_samples > 0) {
        sampled_dataset = rbind(matrix(0, wild_type_samples,11), sampled_dataset)
    }

    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }

    return(sampled_dataset)
}

# sample multiple biopses from the polyclonal-high tree model
sample.multiple.biopses.polyclonal.high <- function (samples_num,
    clones_probabilities,
    nodes_probabilities = NA,
    e_pos,
    e_neg,
    wild_type) {
    
    wild_type_samples = round(samples_num * wild_type)
    if (wild_type > 0 && wild_type_samples == 0) {
        wild_type_samples = 1
    }
    
    res.matrix = matrix(NA, ncol = 17, nrow = (samples_num - wild_type_samples))
    sampled_dataset = sapply(res.matrix[,1], function(x){
        sample.multiple.biospes.from.single.cells(clones_probabilities,
            "high",
            nodes_probabilities,
            0,
            0)
    })
    sampled_dataset = t(sampled_dataset)
    
    if(wild_type_samples > 0) {
        sampled_dataset = rbind(matrix(0, wild_type_samples, 17), sampled_dataset)
    }

    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }

    return(sampled_dataset)
}

# sample multiple biopses from a random single cell tree model
sample.random.multiple.biopses <- function (samples_num,
    e_pos,
    e_neg,
    nodes,
    min_significance,
    max_significance,
    samples_significance,
    wild_type,
    true_tree) {
    
    wild_type_samples = round(samples_num * wild_type)
    if (wild_type > 0 && wild_type_samples == 0) {
        wild_type_samples = 1
    }
    
    # generate a random single cell tree
    random_tree = sample.random.single.cells.significant.tree(nodes,min_significance,max_significance,samples_significance,true_tree)
    
    # as an heuristics, one bulk sample is the mix of nodes/2 single cells samples
    clones_per_sample = round(nodes/2)
    
    res.matrix = matrix(NA, ncol=nodes, nrow=(samples_num - wild_type_samples))
    res = sapply(res.matrix[,1], function(x) {
    	sample.random.single.cells.given.tree(clones_per_sample,
    	                            nodes,
    	                            0,
    	                            0,
    	                            random_tree)
    })
    
    sampled_dataset = t(res)
    
    if (wild_type_samples > 0) {
        sampled_dataset = rbind(matrix(0, wild_type_samples,nodes), sampled_dataset)
    }

    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }
    
    # save the results
    res = list( random_tree = random_tree, sampled_dataset = sampled_dataset )
    
    return(res)
    
}

# apply e_pos and e_neg noise to a sample
apply.noise.to.sample <- function (my_sample,
    e_pos,
    e_neg) {
    
    # apply the noise entry by entry
    for (i in 1:length(my_sample)) {
        # e_pos is the probability of flipping 0 to 1
        if (my_sample[i] == 0) {
            if (runif(1) < e_pos) {
                my_sample[i] = 1
            }
        } else if (my_sample[i] == 1) {
        # e_neg is the probability of flipping 1 to 0
            if (runif(1) < e_neg) {
                my_sample[i] = 0
            }
        }
    }
    return(my_sample)
}

sample.random.single.cells.forest <- function(  samples_num, e_pos, e_neg, nodes, min_significance = 0.60, max_significance = 0.9, sample_significance = 0.005, true_tree = NULL ) {

    num.trees = sample(c(2,3), 1)
    prob.trees = runif(num.trees)
    prob.trees = prob.trees / sum(prob.trees)
    num.nodes = rep(2, num.trees)

    steps = sort(sample(0:(nodes - (num.trees * 2)), num.trees - 1,replace = FALSE))
    steps = c(steps, nodes - (num.trees * 2))

    prec = 0
    for (i in 1:num.trees) {
        nodes.to.add = steps[i] - prec
        prec = steps[i]
        num.nodes[i] = num.nodes[i] + nodes.to.add
    }

    first.node = 0
    last.node = 0

    empty.graph = matrix(0, nodes, nodes)
    colnames(empty.graph) = paste0('node_', 1:nodes)
    rownames(empty.graph) =paste0('node_', 1:nodes)
    final.true.tree = graph.adjacency(empty.graph)
    final.probs = matrix(0, 0, 1)
    final.values = matrix(0, 0, 0)

    for (i in 1:num.trees) {
        tree = sample.random.single.cells(samples_num, 0, 0, num.nodes[[i]], min_significance,max_significance,sample_significance)
        samples = tree$random_tree$dataset_samples
        true.tree = tree$random_tree$structure
        first.node = last.node + 1
        last.node = first.node + num.nodes[[i]] - 1
        node.names = paste0('node_', first.node:last.node)
        colnames(true.tree) = node.names
        rownames(true.tree) = node.names
        colnames(samples) = c(node.names, 'Probs')
        #print(samples)
        #print(prob.trees[i])
        samples[, 'Probs'] = samples[, 'Probs'] * prob.trees[i]
        probs = samples[, 'Probs', drop = FALSE]
        #print(probs)
        final.probs = rbind(final.probs, probs)
        values = samples[, -ncol(samples), drop = FALSE]

        left = rbind(final.values, matrix(0, ncol = ncol(final.values), nrow = nrow(values)))
        right = rbind(matrix(0, ncol = ncol(values), nrow = nrow(final.values)), values)
        final.values = cbind(left, right)

        #print(values)
        #print(samples)
        #print(true.tree)
        true.tree = graph.adjacency(true.tree)
        final.true.tree = graph.union(final.true.tree, true.tree)
    }
    #plot(final.true.tree)
    #print(final.probs)
    #print(final.values)
    final.true.tree = get.adjacency(final.true.tree, sparse = FALSE)
    
    rownames(final.probs) = paste0('sample_', 1:nrow(final.probs))
    rownames(final.values) = paste0('sample_', 1:nrow(final.values))
    final.dataset.samples = cbind(final.values, final.probs)
    final.dataset = final.values[sample(1:nrow(final.values),
        size = samples_num,
        replace = TRUE,
        prob = final.probs),]
    rownames(final.dataset) = paste0("sample_",1:samples_num)

    # apply the noise
    for(i in 1:nrow(final.dataset)) {
        final.dataset[i,] = apply.noise.to.sample(final.dataset[i,], e_pos, e_neg)
    }

    random_tree = list(structure = final.true.tree, dataset_samples = final.dataset.samples)
    res = list( random_tree = random_tree, sampled_dataset = final.dataset )
    return(res)

}



# sample multiple biopses from a random single cell tree model
sample.random.multiple.biopses.forest <- function (samples_num,
    e_pos,
    e_neg,
    nodes,
    min_significance,
    max_significance,
    samples_significance,
    wild_type) {
    
    wild_type_samples = round(samples_num * wild_type)
    if (wild_type > 0 && wild_type_samples == 0) {
        wild_type_samples = 1
    }
    
    # generate a random single cell tree
    random_tree = sample.random.single.cells.forest(samples_num,
        e_pos,
        e_neg,
        nodes,
        min_significance,
        max_significance,
        samples_significance)$random_tree
    
    # as an heuristics, one bulk sample is the mix of nodes/2 single cells samples
    clones_per_sample = round(nodes/2)
    
    res.matrix = matrix(NA, ncol=nodes, nrow=(samples_num - wild_type_samples))
    res = sapply(res.matrix[,1], function(x) {
        sample.random.single.cells.given.tree(clones_per_sample,
                                    nodes,
                                    0,
                                    0,
                                    random_tree)
    })
    
    sampled_dataset = t(res)
    
    if (wild_type_samples > 0) {
        sampled_dataset = rbind(matrix(0, wild_type_samples,nodes), sampled_dataset)
    }

    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }
    
    # save the results
    res = list( random_tree = random_tree, sampled_dataset = sampled_dataset )
    
    return(res)
    
}


sample.random.single.cells.forest.fixed.tree <- function(  samples_num, e_pos, e_neg, min_significance = 0.60, max_significance = 0.9, sample_significance = 0.005) {

    nodes = 7
    num.trees = 2
    num.nodes = c(4,3)
    prob.trees = c(0.6, 0.4)

    first.tree = matrix(0, 4, 4)
    first.tree[1,2] = 1
    first.tree[1,3] = 1
    first.tree[3,4] = 1
    second.tree = matrix(0, 3, 3)
    second.tree[1,2] = 1
    second.tree[2,3] = 1

    true_tree = list(first.tree, second.tree)

    first.node = 0
    last.node = 0

    empty.graph = matrix(0, nodes, nodes)
    colnames(empty.graph) = paste0('node_', 1:nodes)
    rownames(empty.graph) =paste0('node_', 1:nodes)
    final.true.tree = graph.adjacency(empty.graph)
    final.probs = matrix(0, 0, 1)
    final.values = matrix(0, 0, 0)

    for (i in 1:num.trees) {
        tree = sample.random.single.cells(samples_num, 0, 0, num.nodes[[i]], min_significance,max_significance,sample_significance, true_tree[[i]])
        samples = tree$random_tree$dataset_samples
        true.tree = tree$random_tree$structure
        first.node = last.node + 1
        last.node = first.node + num.nodes[[i]] - 1
        node.names = paste0('node_', first.node:last.node)
        colnames(true.tree) = node.names
        rownames(true.tree) = node.names
        colnames(samples) = c(node.names, 'Probs')
        #print(samples)
        #print(prob.trees[i])
        samples[, 'Probs'] = samples[, 'Probs'] * prob.trees[i]
        probs = samples[, 'Probs', drop = FALSE]
        #print(probs)
        final.probs = rbind(final.probs, probs)
        values = samples[, -ncol(samples), drop = FALSE]

        left = rbind(final.values, matrix(0, ncol = ncol(final.values), nrow = nrow(values)))
        right = rbind(matrix(0, ncol = ncol(values), nrow = nrow(final.values)), values)
        final.values = cbind(left, right)

        #print(values)
        #print(samples)
        #print(true.tree)
        true.tree = graph.adjacency(true.tree)
        final.true.tree = graph.union(final.true.tree, true.tree)
    }
    #plot(final.true.tree)
    #print(final.probs)
    #print(final.values)
    final.true.tree = get.adjacency(final.true.tree, sparse = FALSE)
    
    rownames(final.probs) = paste0('sample_', 1:nrow(final.probs))
    rownames(final.values) = paste0('sample_', 1:nrow(final.values))
    final.dataset.samples = cbind(final.values, final.probs)
    final.dataset = final.values[sample(1:nrow(final.values),
        size = samples_num,
        replace = TRUE,
        prob = final.probs),]
    rownames(final.dataset) = paste0("sample_",1:samples_num)

    # apply the noise
    for(i in 1:nrow(final.dataset)) {
        final.dataset[i,] = apply.noise.to.sample(final.dataset[i,], e_pos, e_neg)
    }

    random_tree = list(structure = final.true.tree, dataset_samples = final.dataset.samples)
    res = list( random_tree = random_tree, sampled_dataset = final.dataset )
    return(res)

}

# sample multiple biopses from a random single cell tree model
sample.random.multiple.biopses.forest.fixed.tree <- function (samples_num,
    e_pos,
    e_neg,
    min_significance,
    max_significance,
    samples_significance,
    wild_type) {
    
    wild_type_samples = round(samples_num * wild_type)
    if (wild_type > 0 && wild_type_samples == 0) {
        wild_type_samples = 1
    }
    
    # generate a random single cell tree
    random_tree = sample.random.single.cells.forest.fixed.tree(samples_num,
        e_pos,
        e_neg,
        min_significance,
        max_significance,
        samples_significance)$random_tree
    
    # as an heuristics, one bulk sample is the mix of nodes/2 single cells samples
    nodes = 7
    clones_per_sample = round(nodes/2)
    
    res.matrix = matrix(NA, ncol=nodes, nrow=(samples_num - wild_type_samples))
    res = sapply(res.matrix[,1], function(x) {
        sample.random.single.cells.given.tree(clones_per_sample,
                                    nodes,
                                    0,
                                    0,
                                    random_tree)
    })
    
    sampled_dataset = t(res)
    
    if (wild_type_samples > 0) {
        sampled_dataset = rbind(matrix(0, wild_type_samples,nodes), sampled_dataset)
    }

    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }
    
    # save the results
    res = list( random_tree = random_tree, sampled_dataset = sampled_dataset )
    
    return(res)
    
}



# generate a dataset of random single cell trees
sample.random.single.cells.convergent <- function ( samples_num, e_pos, e_neg, nodes, min_significance = 0.60, max_significance = 0.9, sample_significance = 0.005, true_tree = NULL ) {
    # generate a random tree with minimum significance per samples
    random_tree = NULL
    while(is.null(random_tree) || 
            min(random_tree$dataset_samples[,"Probs"])<=sample_significance || 
            max(random_tree$dataset_samples[,"Probs"])>=(1-sample_significance)) {
                

        # keep generating random trees until I get a valid one
        random_tree = generate.random.single.cell.dataset.convergent(nodes,min_significance,max_significance,true_tree)
        
    }
    
    dataset = random_tree$dataset_samples[,1:nodes]
    samples_probabilities = random_tree$dataset_samples[,nodes+1]
    
    # sample the dataset
    sampled_dataset = dataset[sample(1:nrow(dataset),
        size = samples_num,
        replace = TRUE,
        prob = samples_probabilities),]
    
    # apply the noise
    for(i in 1:nrow(sampled_dataset)) {
        sampled_dataset[i,] = apply.noise.to.sample(sampled_dataset[i,], e_pos, e_neg)
    }
    rownames(sampled_dataset) = paste0("sample_",1:samples_num)
    
    # save the results
    res = list( random_tree = random_tree, sampled_dataset = sampled_dataset )
    
    return(res)
    
}
