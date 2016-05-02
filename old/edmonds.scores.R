perform.scores.edmonds = function( tronco.object, score.type ) {
    #print(score.type)
    dataset = as.genotypes(tronco.object)
    adj.matrix = tronco.object$adj.matrix.prima.facie
    marginal.probs = as.marginal.probs(tronco.object)[[1]]
    joint.probs = as.joint.probs(tronco.object)[[1]]

    adj.matrix.prima.facie = adj.matrix
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
       
    # set at most one parent per node based on mutual information
    for (i in 1:ncol(adj.matrix)) {
            
        # consider the parents of i
        curr_parents = which(adj.matrix[,i] == 1)
        
        # if I have more then one valid parent
        if (length(curr_parents) > 1) {
            
            # find the best parent
            curr_best_parent = -1
            curr_best_score = -1
            for (j in curr_parents) {
                # mutual information with bootstrap
                if(score.type=="MI") {
                    #print(joint.probs[i,j])
                    #print(marginal.probs[i])
                    #print(marginal.probs[j])
                    if (joint.probs[i,j] > 0) {
                        p_j_i = joint.probs[i,j] * log( joint.probs[i,j] / (marginal.probs[i] * marginal.probs[j]))
                    } else {
                        p_j_i = 1 # marginal.probs[i] * log( marginal.probs[i] / (marginal.probs[i] * marginal.probs[j]))
                    }
                    #print(p_j_i)
                    if ((marginal.probs[i] - joint.probs[i,j]) > 0 && joint.probs[i,j] != -1) {
                        p_not_j_i = (marginal.probs[i] - joint.probs[i,j]) * log( (marginal.probs[i] - joint.probs[i,j]) / (marginal.probs[i] * (1 - marginal.probs[j])))
                    } else if ((marginal.probs[i] - joint.probs[i,j]) == 0 && joint.probs[i,j] != -1) {
                        p_not_j_i = 0
                    } else {
                        p_not_j_i = 0
                    }
                    #print(p_not_j_i)
                    if ((marginal.probs[j] - joint.probs[i,j]) > 0 && joint.probs[i,j] != -1) {
                        p_not_i_j = (marginal.probs[j] - joint.probs[i,j]) * log( (marginal.probs[j] - joint.probs[i,j]) / ((1 - marginal.probs[i]) * marginal.probs[j]))
                    } else if ((marginal.probs[j] - joint.probs[i,j]) == 0 && joint.probs[i,j] != -1) {
                        p_not_i_j = 0
                    } else {
                        p_not_i_j = 0
                    }
                    #print(p_not_i_j)
                    if (joint.probs[i,j] > 0) {
                        p_not_i_not_j = (1 - marginal.probs[i] - marginal.probs[j] + joint.probs[i,j]) * log((1 - marginal.probs[i] - marginal.probs[j] + joint.probs[i,j]) / ((1 - marginal.probs[i]) * (1 - marginal.probs[j])))
                    } else {
                        p_not_i_not_j = 0 # (1 - marginal.probs[i] - marginal.probs[j] + marginal.probs[i]) * log((1 - marginal.probs[i] - marginal.probs[j] + marginal.probs[i]) / ((1 - marginal.probs[i]) * (1 - marginal.probs[j])))
                    }
                    #print(p_not_i_not_j)
                    new_score = p_j_i + p_not_i_j + p_not_j_i + p_not_i_not_j
                }
                # causal mutual information
                else if(score.type=="CMI") {

                    if (joint.probs[i,j] > 0) {
                        p_j_i = joint.probs[i,j] * log( joint.probs[i,j] / (marginal.probs[i] * marginal.probs[j]))
                    } else {
                        p_j_i = 1 # marginal.probs[i] * log( marginal.probs[i] / (marginal.probs[i] * marginal.probs[j]))
                    }

                    if ((marginal.probs[j] - joint.probs[i,j]) > 0 && joint.probs[i,j] != -1) {
                        p_not_i_j = (marginal.probs[j] - joint.probs[i,j]) * log( (marginal.probs[j] - joint.probs[i,j]) / ((1 - marginal.probs[i]) * marginal.probs[j]))
                    } else if ((marginal.probs[j] - joint.probs[i,j]) == 0 && joint.probs[i,j] != -1) {
                        p_not_i_j = 0
                    } else {
                        p_not_i_j = 0
                    }

                    new_score = p_j_i + p_not_i_j
                }
                # causal mutual information - positive and neutral
                else if(score.type=="CMI2") {

                    if (joint.probs[i,j] > 0) {
                        p_j_i = joint.probs[i,j] * log( joint.probs[i,j] / (marginal.probs[i] * marginal.probs[j]))
                    } else {
                        p_j_i = 1 # marginal.probs[i] * log( marginal.probs[i] / (marginal.probs[i] * marginal.probs[j]))
                    }

                    if (joint.probs[i,j] > 0) {
                        p_not_i_not_j = (1 - marginal.probs[i] - marginal.probs[j] + joint.probs[i,j]) * log((1 - marginal.probs[i] - marginal.probs[j] + joint.probs[i,j]) / ((1 - marginal.probs[i]) * (1 - marginal.probs[j])))
                    } else {
                        p_not_i_not_j = 0 #(1 - marginal.probs[i] - marginal.probs[j] + marginal.probs[i]) * log((1 - marginal.probs[i] - marginal.probs[j] + marginal.probs[i]) / ((1 - marginal.probs[i]) * (1 - marginal.probs[j])))
                    }

                    new_score = p_j_i + p_not_i_not_j
                }
                # causal mutual information - positive, neutral and no parent
                else if(score.type=="CMI3") {

                    if (joint.probs[i,j] > 0) {
                        p_j_i = joint.probs[i,j] * log( joint.probs[i,j] / (marginal.probs[i] * marginal.probs[j]))
                    } else {
                        p_j_i = 1 #marginal.probs[i] * log( marginal.probs[i] / (marginal.probs[i] * marginal.probs[j]))
                    }

                    if ((marginal.probs[j] - joint.probs[i,j]) > 0 && joint.probs[i,j] != -1) {
                        p_not_i_j = (marginal.probs[j] - joint.probs[i,j]) * log( (marginal.probs[j] - joint.probs[i,j]) / ((1 - marginal.probs[i]) * marginal.probs[j]))
                    } else if ((marginal.probs[j] - joint.probs[i,j]) == 0 && joint.probs[i,j] != -1) {
                        p_not_i_j = 0
                    } else {
                        p_not_i_j = 0
                    }

                    if (joint.probs[i,j] > 0) {
                        p_not_i_not_j = (1 - marginal.probs[i] - marginal.probs[j] + joint.probs[i,j]) * log((1 - marginal.probs[i] - marginal.probs[j] + joint.probs[i,j]) / ((1 - marginal.probs[i]) * (1 - marginal.probs[j])))
                    } else {
                        p_not_i_not_j = 0 #(1 - marginal.probs[i] - marginal.probs[j] + marginal.probs[i]) * log((1 - marginal.probs[i] - marginal.probs[j] + marginal.probs[i]) / ((1 - marginal.probs[i]) * (1 - marginal.probs[j])))
                    }

                    new_score = p_j_i + p_not_i_not_j + p_not_i_j
                }
                # pointwise mutual information
                else if(score.type=="OR") {

                    if (joint.probs[i,j] > 0) {
                        p_j_i = log( joint.probs[i,j] / (marginal.probs[i] * marginal.probs[j]))
                    } else {
                        p_j_i = Inf #log( marginal.probs[i] / (marginal.probs[i] * marginal.probs[j]))
                    }

                    new_score = p_j_i
                }
                else if(score.type=="OR2") {

                    if (joint.probs[i,j] > 0) {
                        p_j_i = log( joint.probs[i,j] / (marginal.probs[i] * marginal.probs[j]))

                        p_causal_1 = joint.probs[i,j] / marginal.probs[j]
                        p_causal_2 = (marginal.probs[i] - joint.probs[i,j]) / (1 - marginal.probs[j])
                        p_causal = (p_causal_1 - p_causal_2) / (p_causal_1 + p_causal_2)
                        p_j_i = p_j_i * p_causal

                    } else {
                        p_j_i = Inf #log( marginal.probs[i] / (marginal.probs[i] * marginal.probs[j]))
                    }

                    new_score = p_j_i
                }

                # correlation pearson
                #else if(score.type=="pearson") {
                #    scores = cor(dataset,method="pearson")
                #    new_score = scores[i,j]
                #    if(is.na(new_score)) {
                #        new_score = -1
                #    }
                #}
                ## correlation kendall
                #else if(score.type=="kendall") {
                #    scores = cor(dataset,method="kendall")
                #    new_score = scores[i,j]
                #    if(is.na(new_score)) {
                #        new_score = -1
                #    }
                #}
                ## correlation spearman
                #else if(score.type=="spearman") {
                #    scores = cor(dataset,method="spearman")
                #    new_score = scores[i,j]
                #    if(is.na(new_score)) {
                #        new_score = -1
                #    }
                #}
                # causal pointwise mutual information 1
                #else if(score.type=="CAUSAL_OR_1") {
                #    
                #    p_i_given_j = joint.probs[i,j]/marginal.probs[j]
                #    p_i_given_not_j = (marginal.probs[i]-joint.probs[i,j])/(1-marginal.probs[j])
                #    causal_score = (p_i_given_j - p_i_given_not_j) / (p_i_given_j + p_i_given_not_j)
                #    
                #    correlation_score = log(joint.probs[i,j]/(marginal.probs[i]*marginal.probs[j]))
                #    
                #    new_score = causal_score * correlation_score
                #    
                #}
                ## causal pointwise mutual information 2
                #else if(score.type=="CAUSAL_OR_2") {
                #    
                #    p_i_given_j = joint.probs[i,j]/marginal.probs[j]
                #    p_i_given_not_j = (marginal.probs[i]-joint.probs[i,j])/(1-marginal.probs[j])
                #    causal_score = (p_i_given_j - p_i_given_not_j) / (p_i_given_j + p_i_given_not_j)
                #    
                #    correlation_score_2 = joint.probs[i,j]/(marginal.probs[i]*marginal.probs[j])
                #    
                #    new_score = causal_score * correlation_score_2
                #    
                #}
                ## causal pointwise mutual information 3
                #else if(score.type=="CAUSAL_OR_3") {
                #    
                #    p_i_given_j = joint.probs[i,j]/marginal.probs[j]
                #    p_i_given_not_j = (marginal.probs[i]-joint.probs[i,j])/(1-marginal.probs[j])
                #    
                #    new_score = (p_i_given_j - p_i_given_not_j + joint.probs[i,j]) / (p_i_given_j + p_i_given_not_j + marginal.probs[i]*marginal.probs[j])
                #    
                #} # causal pointwise mutual information
                #else if(score.type=="CAUSAL_OR_4") {
                #    
                #    # pointwise mutual information for (1,1) as positive case to measure positive dependence
                #    score_1 = log(joint.probs[i,j]/(marginal.probs[i]*marginal.probs[j]))
                #    
                #    # negated normalized pointwise mutual information for (0,1) as negative case to measure the probability raising
                #    score_2 = -log((marginal.probs[i]-joint.probs[i,j])/(marginal.probs[i]*(1-marginal.probs[j])))/(-log((marginal.probs[i]-joint.probs[i,j])))
                #    
                #    if (is.nan(score_2)) {
                #        score_2 = 1
                #    }
                #    new_score = score_1 * score_2
                #    
                #}
                else {
                    stop('error!')
                }

                if (new_score > curr_best_score) {
                    curr_best_parent = j
                    curr_best_score = new_score
                }
            
            }
            
            # set the best parent
            for (j in curr_parents) {
                if (j != curr_best_parent) {
                    adj.matrix[j,i] = 0
                }
            }
        }
    }
    adj.matrix.fit = adj.matrix
    
    ## Save the results and return them.
    
    return(adj.matrix.fit)

}
