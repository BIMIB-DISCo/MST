
# perform the reconstructions
run.reconstructions <- function( dataset, true_tree ) {
    
    results = NULL
    
    # create a TRONCO object of the dataset
    data = import.genotypes(dataset)
    
    ## performs the reconstructions with CAPRI loglik, aic and bic
    res = NULL
    #res = tronco.capri(data,regularization=c("loglik","aic","bic"), silent = TRUE)
    res = tronco.capri(data,regularization=c("aic","bic"), silent = TRUE)
    #adj.matrix.capri.loglik = as.adj.matrix(res,model="capri_loglik")
    #results.capri.loglik = getStats(true_tree,adj.matrix.capri.loglik[["capri_loglik"]])
    adj.matrix.capri.aic = as.adj.matrix(res,model="capri_aic")
    results.capri.aic = getStats(true_tree,adj.matrix.capri.aic[["capri_aic"]])
    adj.matrix.capri.bic = as.adj.matrix(res,model="capri_bic")
    results.capri.bic = getStats(true_tree,adj.matrix.capri.bic[["capri_bic"]])
    #capri = list(loglik.adj=adj.matrix.capri.loglik,loglik.res=results.capri.loglik,aic.adj=adj.matrix.capri.aic,aic.res=results.capri.aic,bic.adj=adj.matrix.capri.bic,bic.res=results.capri.bic)
    capri = list(aic.adj=adj.matrix.capri.aic,aic.res=results.capri.aic,bic.adj=adj.matrix.capri.bic,bic.res=results.capri.bic)
    results[["capri"]] = capri
    
    # performs the reconstructions with CAPRESE
    res = NULL
    res = tronco.caprese(data, silent = TRUE)
    adj.matrix.caprese = as.adj.matrix(res,model="caprese")
    results.caprese = getStats(true_tree,adj.matrix.caprese[["caprese"]])
    caprese = list(adj.caprese=adj.matrix.caprese,caprese.res=results.caprese)
    results[["caprese"]] = caprese
    
    # performs the reconstructions with Edmonds no_reg, loglik, aic and bic
    res = NULL
    #res = tronco.mst.edmonds(data,regularization=c("no_reg","loglik","aic","bic"))
    res = tronco.mst.edmonds(data,regularization=c("no_reg"))
    adj.matrix.edmonds.no.reg = as.adj.matrix(res,model="edmonds_no_reg")
    results.edmonds.no.reg = getStats(true_tree,adj.matrix.edmonds.no.reg[["edmonds_no_reg"]])  
    #adj.matrix.edmonds.loglik = as.adj.matrix(res,model="edmonds_loglik")
    #results.edmonds.loglik = getStats(true_tree,adj.matrix.edmonds.loglik[["edmonds_loglik"]])
    #adj.matrix.edmonds.aic = as.adj.matrix(res,model="edmonds_aic")
    #results.edmonds.aic = getStats(true_tree,adj.matrix.edmonds.aic[["edmonds_aic"]])
    #adj.matrix.edmonds.bic = as.adj.matrix(res,model="edmonds_bic")
    #results.edmonds.bic = getStats(true_tree,adj.matrix.edmonds.bic[["edmonds_bic"]])
    #edmonds = list(no.reg.adj=adj.matrix.edmonds.no.reg,no.reg.res=results.edmonds.no.reg,loglik.adj=adj.matrix.edmonds.loglik,loglik.res=results.edmonds.loglik,aic.adj=adj.matrix.edmonds.aic,aic.res=results.edmonds.aic,bic.adj=adj.matrix.edmonds.bic,bic.res=results.edmonds.bic)
    edmonds = list(no.reg.adj=adj.matrix.edmonds.no.reg,no.reg.res=results.edmonds.no.reg)
    results[["edmonds"]] = edmonds
    
    # performs the reconstructions with Chow Liu loglik, aic and bic
    res = NULL
    #res = tronco.mst.chowliu(data,regularization=c("loglik","aic","bic"))   
    res = tronco.mst.chowliu(data,regularization=c("loglik"))   
    adj.matrix.chowliu.loglik = as.adj.matrix(res,model="chow_liu_loglik")
    results.chowliu.loglik = getStats(true_tree,adj.matrix.chowliu.loglik[["chow_liu_loglik"]])
    #adj.matrix.chowliu.aic = as.adj.matrix(res,model="chow_liu_aic")
    #results.chowliu.aic = getStats(true_tree,adj.matrix.chowliu.aic[["chow_liu_aic"]])
    #adj.matrix.chowliu.bic = as.adj.matrix(res,model="chow_liu_bic")
    #results.chowliu.bic = getStats(true_tree,adj.matrix.chowliu.bic[["chow_liu_bic"]])
    chowliu = list(loglik.adj=adj.matrix.chowliu.loglik,loglik.res=results.chowliu.loglik)
    results[["chowliu"]] = chowliu
    
    # performs the reconstructions with Prim no_reg, loglik, aic and bic
    res = NULL
    #res = tronco.mst.prim(data,regularization=c("no_reg","loglik","aic","bic"))
    res = tronco.mst.prim(data,regularization=c("no_reg"))
    adj.matrix.prim.no.reg = as.adj.matrix(res,model="prim_no_reg")
    results.prim.no.reg = getStats(true_tree,adj.matrix.prim.no.reg[["prim_no_reg"]])
    #adj.matrix.prim.loglik = as.adj.matrix(res,model="prim_loglik")
    #results.prim.loglik = getStats(true_tree,adj.matrix.prim.loglik[["prim_loglik"]])
    #adj.matrix.prim.aic = as.adj.matrix(res,model="prim_aic")
    #results.prim.aic = getStats(true_tree,adj.matrix.prim.aic[["prim_aic"]])
    #adj.matrix.prim.bic = as.adj.matrix(res,model="prim_bic")
    #results.prim.bic = getStats(true_tree,adj.matrix.prim.bic[["prim_bic"]])
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



expand.input <- function(datasets, true_tree, seed, cores) {
    cat('Using', cores, 'cores via "parallel" \n')
    cl = makeCluster(cores, outfile='')
    clusterEvalQ(cl, library(TRONCO))
    clusterExport(cl, 'run.reconstructions')
    clusterExport(cl, 'getStats')
    clusterSetRNGStream(cl, iseed = seed)

    for(i in 1:nrow(datasets)) {
        for (j in 1:ncol(datasets)) {
            cat((((i - 1) * ncol(datasets)) + j) , '/', nrow(datasets) * ncol(datasets), '\n')
            single.experiment = datasets[[i,j]]       
            results = parLapply(cl, single.experiment, function(x, true_tree){
                run.reconstructions(x$dataset, true_tree)
            }, true_tree)

            for (k in 1:length(results)) {
                datasets[[i,j]][[k]]$reconstructions = results[[k]]
            }
        }
    }
    stopCluster(cl)
    return(datasets)

}

expand.random.input <- function(datasets, seed, cores) {
    cat('Using', cores, 'cores via "parallel" \n')
    cl = makeCluster(cores)
    clusterEvalQ(cl, library(TRONCO))
    clusterExport(cl, 'run.reconstructions')
    clusterExport(cl, 'getStats')
    clusterSetRNGStream(cl, iseed = seed)

    for(i in 1:nrow(datasets)) {
        for (j in 1:ncol(datasets)) {
            cat((((i - 1) * ncol(datasets)) + j) , '/', nrow(datasets) * ncol(datasets), '\n')
            single.experiment = datasets[[i,j]]       
            results = parLapply(cl, single.experiment, function(x){
                run.reconstructions(x$dataset, x$true_tree$structure)
            })

            for (k in 1:length(results)) {
                datasets[[i,j]][[k]]$reconstructions = results[[k]]
            }
        }
    }
    stopCluster(cl)
    return(datasets)

}


