
# perform the reconstructions
run.reconstructions <- function( dataset, true_tree, epos, eneg ) {
    
    results = NULL
    
    # create a TRONCO object of the dataset
    data = import.genotypes(dataset)
    
    ## performs the reconstructions with CAPRI loglik, aic and bic
    res = NULL
    res = tronco.capri(data,regularization=c("loglik","aic","bic"), silent = TRUE, pvalue = 0.1, nboot = 1000, epos = epos, eneg = eneg)
    adj.matrix.capri.loglik = as.adj.matrix(res,model="capri_loglik")
    results.capri.loglik = getStats(true_tree,adj.matrix.capri.loglik[["capri_loglik"]])
    adj.matrix.capri.aic = as.adj.matrix(res,model="capri_aic")
    results.capri.aic = getStats(true_tree,adj.matrix.capri.aic[["capri_aic"]])
    adj.matrix.capri.bic = as.adj.matrix(res,model="capri_bic")
    results.capri.bic = getStats(true_tree,adj.matrix.capri.bic[["capri_bic"]])
    capri = list(loglik.adj=adj.matrix.capri.loglik,loglik.res=results.capri.loglik,
        aic.adj=adj.matrix.capri.aic,aic.res=results.capri.aic,
        bic.adj=adj.matrix.capri.bic,bic.res=results.capri.bic)
    results[["capri"]] = capri
    
    # performs the reconstructions with CAPRESE
    res = NULL
    res = tronco.caprese(data, silent = TRUE, epos = epos, eneg = eneg)
    adj.matrix.caprese = as.adj.matrix(res,model="caprese")
    results.caprese = getStats(true_tree,adj.matrix.caprese[["caprese"]])
    caprese = list(adj.caprese=adj.matrix.caprese,caprese.res=results.caprese)
    results[["caprese"]] = caprese
    
    # performs the reconstructions with Edmonds no_reg, loglik, aic and bic
    res = NULL
    res = tronco.mst.edmonds(data,regularization=c("no_reg"), score=c('entropy', 'pmi', 'cpmi'), pvalue = 0.1, nboot = 1000, epos = epos, eneg = eneg, silent = TRUE)
    adj.matrix.edmonds.entropy.no.reg = as.adj.matrix(res,model="edmonds_no_reg_entropy")
    results.edmonds.entropy.no.reg = getStats(true_tree,adj.matrix.edmonds.entropy.no.reg[["edmonds_no_reg_entropy"]])

    adj.matrix.edmonds.pmi.no.reg = as.adj.matrix(res,model="edmonds_no_reg_pmi")
    results.edmonds.pmi.no.reg = getStats(true_tree,adj.matrix.edmonds.pmi.no.reg[["edmonds_no_reg_pmi"]])

    adj.matrix.edmonds.cpmi.no.reg = as.adj.matrix(res,model="edmonds_no_reg_cpmi")
    results.edmonds.cpmi.no.reg = getStats(true_tree,adj.matrix.edmonds.cpmi.no.reg[["edmonds_no_reg_cpmi"]])

    edmonds = list(entropy.no.reg.adj=adj.matrix.edmonds.entropy.no.reg, entropy.no.reg.res=results.edmonds.entropy.no.reg,
        pmi.no.reg.adj = adj.matrix.edmonds.pmi.no.reg, pmi.no.reg.res = results.edmonds.pmi.no.reg,
        cpmi.no.reg.adj = adj.matrix.edmonds.cpmi.no.reg, cpmi.no.reg.res = results.edmonds.cpmi.no.reg)
    results[["edmonds"]] = edmonds

    # performs the reconstructions with Gabow
    res = NULL
    res = tronco.mst.gabow(data,regularization=c("no_reg"), score=c('entropy', 'pmi', 'cpmi', 'mi'), pvalue = 0.1, nboot = 1000, epos = epos, eneg = eneg, silent = TRUE)
    adj.matrix.gabow.entropy.no.reg = as.adj.matrix(res,model="gabow_no_reg_entropy")
    results.gabow.entropy.no.reg = getStats(true_tree,adj.matrix.gabow.entropy.no.reg[["gabow_no_reg_entropy"]])

    adj.matrix.gabow.pmi.no.reg = as.adj.matrix(res,model="gabow_no_reg_pmi")
    results.gabow.pmi.no.reg = getStats(true_tree,adj.matrix.gabow.pmi.no.reg[["gabow_no_reg_pmi"]])

    adj.matrix.gabow.cpmi.no.reg = as.adj.matrix(res,model="gabow_no_reg_cpmi")
    results.gabow.cpmi.no.reg = getStats(true_tree,adj.matrix.gabow.cpmi.no.reg[["gabow_no_reg_cpmi"]])

    adj.matrix.gabow.mi.no.reg = as.adj.matrix(res,model="gabow_no_reg_mi")
    results.gabow.mi.no.reg = getStats(true_tree,adj.matrix.gabow.mi.no.reg[["gabow_no_reg_mi"]])

    gabow = list(entropy.no.reg.adj=adj.matrix.gabow.entropy.no.reg, entropy.no.reg.res=results.gabow.entropy.no.reg,
        pmi.no.reg.adj = adj.matrix.gabow.pmi.no.reg, pmi.no.reg.res = results.gabow.pmi.no.reg,
        cpmi.no.reg.adj = adj.matrix.gabow.cpmi.no.reg, cpmi.no.reg.res = results.gabow.cpmi.no.reg,
        mi.no.reg.adj = adj.matrix.gabow.mi.no.reg, mi.no.reg.res = results.gabow.mi.no.reg)
    results[["gabow"]] = gabow

    # performs the reconstructions with Gabow no raising
    res = NULL
    res = tronco.mst.gabow(data,regularization=c("no_reg"), score=c('entropy', 'pmi', 'cpmi', 'mi'), pvalue = 0.1, nboot = 1000, epos = epos, eneg = eneg, do.raising = FALSE, silent = TRUE)
    
    adj.matrix.gabow.no.raising.entropy.no.reg = as.adj.matrix(res,model="gabow_no_reg_entropy")
    results.gabow.no.raising.entropy.no.reg = getStats(true_tree,adj.matrix.gabow.no.raising.entropy.no.reg[["gabow_no_reg_entropy"]])

    adj.matrix.gabow.no.raising.pmi.no.reg = as.adj.matrix(res,model="gabow_no_reg_pmi")
    results.gabow.no.raising.pmi.no.reg = getStats(true_tree,adj.matrix.gabow.no.raising.pmi.no.reg[["gabow_no_reg_pmi"]])

    adj.matrix.gabow.no.raising.cpmi.no.reg = as.adj.matrix(res,model="gabow_no_reg_cpmi")
    results.gabow.no.raising.cpmi.no.reg = getStats(true_tree,adj.matrix.gabow.no.raising.cpmi.no.reg[["gabow_no_reg_cpmi"]])

    adj.matrix.gabow.no.raising.mi.no.reg = as.adj.matrix(res,model="gabow_no_reg_mi")
    results.gabow.no.raising.mi.no.reg = getStats(true_tree,adj.matrix.gabow.no.raising.mi.no.reg[["gabow_no_reg_mi"]])

    gabow = list(no.raising.entropy.no.reg.adj=adj.matrix.gabow.no.raising.entropy.no.reg, no.raising.entropy.no.reg.res=results.gabow.no.raising.entropy.no.reg,
        no.raising.pmi.no.reg.adj = adj.matrix.gabow.no.raising.pmi.no.reg, no.raising.pmi.no.reg.res = results.gabow.no.raising.pmi.no.reg,
        no.raising.cpmi.no.reg.adj = adj.matrix.gabow.no.raising.cpmi.no.reg, no.raising.cpmi.no.reg.res = results.gabow.no.raising.cpmi.no.reg,
        no.raising.mi.no.reg.adj = adj.matrix.gabow.no.raising.mi.no.reg, no.raising.mi.no.reg.res = results.gabow.no.raising.mi.no.reg)
    results[["gabow_no_rising"]] = gabow
      
    # performs the reconstructions with Chow Liu loglik, aic and bic
    res = NULL
    res = tronco.mst.chowliu(data,regularization=c("loglik"), pvalue = 0.1, nboot = 1000, epos = epos, eneg = eneg, silent = TRUE)   
    adj.matrix.chowliu.loglik = as.adj.matrix(res,model="chow_liu_loglik")
    results.chowliu.loglik = getStats(true_tree,adj.matrix.chowliu.loglik[["chow_liu_loglik"]])
    chowliu = list(loglik.adj=adj.matrix.chowliu.loglik,loglik.res=results.chowliu.loglik)
    results[["chowliu"]] = chowliu
    
    # performs the reconstructions with Prim no_reg, loglik, aic and bic
    res = NULL
    res = tronco.mst.prim(data,regularization=c("no_reg"), pvalue = 0.1, nboot = 1000, epos = epos, eneg = eneg, silent = TRUE)
    adj.matrix.prim.no.reg = as.adj.matrix(res,model="prim_no_reg")
    results.prim.no.reg = getStats(true_tree,adj.matrix.prim.no.reg[["prim_no_reg"]])
    prim = list(no.reg.adj=adj.matrix.prim.no.reg,no.reg.res=results.prim.no.reg)
    results[["prim"]] = prim
    
    return(results)
    
}

# function to compute the results at each step
getStats <- function(true_matrix,
    inferred_matrix){
    
    #print(true_matrix)
    #print(inferred_matrix)

    # compute the statistics
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    for (i in 1:nrow(inferred_matrix)) {
        for (j in i:ncol(inferred_matrix)) {
            if (i != j) {
                if (true_matrix[i,j] == 0 && inferred_matrix[i,j] == 0) {
                    tn = tn + 1
                } else if (true_matrix[i,j] == 0 && inferred_matrix[i,j] == 1) {
                    fp = fp + 1
                } else if (true_matrix[i,j] == 1 && inferred_matrix[i,j] == 0) {
                    fn = fn + 1
                } else if (true_matrix[i,j] == 1 && inferred_matrix[i,j] == 1) {
                    tp = tp + 1
                }
            }
        }
    }
    
    # compute the statistics
    accuracy = (tp+tn)/(tp+tn+fp+fn)
    sensitivity = (tp)/(tp+fn)
    specificity = (tn)/(fp+tn)
    hamming_distance = fp + fn
    
    #print(accuray)
    #print(sensitivity)
    #print(specificity)
    
    # check for NAs generated by 0/0 computations
    # this can happen for sensitivity or specificity
    if(is.na(sensitivity)) {
    	sensitivity = 0
    }
    if(is.na(specificity)) {
    	specificity = 0
    }

    # return the results
    results_values = list(accuracy=accuracy,sensitivity=sensitivity,specificity=specificity,hamming_distance=hamming_distance)
    return(results_values)
    
}



expand.input <- function(datasets, true_tree, seed, cores, epos.list, eneg.list) {
    cat('Using', cores, 'cores via "parallel" \n')
    #cl = makeCluster(cores, outfile='')
    #clusterEvalQ(cl, library(TRONCO))
    #clusterExport(cl, 'run.reconstructions')
    #clusterExport(cl, 'getStats')
    #clusterSetRNGStream(cl, iseed = seed)

    for(i in 1:nrow(datasets)) {
        for (j in 1:ncol(datasets)) {
            cat((((i - 1) * ncol(datasets)) + j) , '/', nrow(datasets) * ncol(datasets), '\n')
            for (k in 1:length(datasets[[i,j]])) {
                datasets[[i,j]][[k]]$true_tree$epos = epos.list[[k]]
                datasets[[i,j]][[k]]$true_tree$eneg = eneg.list[[k]]
            }
            single.experiment = datasets[[i,j]]   
            results = lapply(single.experiment, function(x){
                run.reconstructions(x$dataset,
                    true_tree,
                    x$true_tree$epos,
                    x$true_tree$eneg)
            })
            for (k in 1:length(results)) {
                datasets[[i,j]][[k]]$reconstructions = results[[k]]
            }
        }
    }
    #stopCluster(cl)
    return(datasets)

}

