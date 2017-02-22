# perform the inference
load("RData/dataset.missing.data.RData")
library(TRONCO)
library(igraph)
library(bnlearn)

# set the seed
set.seed(744911)

run.reconstructions <- function( dataset, epos, eneg, debug = FALSE, seed = NA, pass.error.rates = TRUE) {

    if (!pass.error.rates) {
        epos = 0
        eneg = 0
    }
    
    results = NULL
    
    silent = !debug

    if (is.na(seed)) {
        seed = as.integer(runif(1) * 10000)
    }

    # create a TRONCO object of the dataset
    data = import.genotypes(dataset)

    # create the igraph structure
    net = empty.graph(colnames(as.genotypes(data)), num = 6)
    categoric.dataset = data.frame(apply(as.genotypes(data), 2, factor))
    for (name in colnames(categoric.dataset)) {
        levels(categoric.dataset[[name]]) = c(0,1)
    }
   
    # performs the reconstructions with CAPRI loglik, aic and bic
    res = NULL
    res = tronco.capri(data,regularization=c("bic"), epos = epos, eneg = eneg, silent = silent, boot.seed = seed)
    adj.matrix.capri.bic = as.adj.matrix(res,model="capri_bic")[['capri_bic']]
    amat(net[[1]]) = adj.matrix.capri.bic
    score = bnlearn::score(net[[1]], categoric.dataset, type='loglik') 
    capri = list(bic.adj=adj.matrix.capri.bic, score = score, model=res)
    results[["capri"]] = capri
    
    # performs the reconstructions with CAPRESE
    res = NULL
    res = tronco.caprese(data, epos = epos, eneg = eneg, silent = silent)
    adj.matrix.caprese = as.adj.matrix(res,model="caprese")[['caprese']]
    amat(net[[2]]) = adj.matrix.caprese
    score = bnlearn::score(net[[2]], categoric.dataset, type='loglik') 
    caprese = list(adj.caprese=adj.matrix.caprese, score = score, model=res)
    results[["caprese"]] = caprese
    
    # performs the reconstructions with Edmonds no_reg, loglik, aic and bic
    res = NULL
    res = tronco.edmonds(data,regularization=c("no_reg"), score=c('pmi'), epos = epos, eneg = eneg, silent = silent, boot.seed = seed)
    adj.matrix.edmonds.pmi.no.reg = as.adj.matrix(res,model="edmonds_no_reg_pmi")[['edmonds_no_reg_pmi']]
    amat(net[[3]]) = adj.matrix.edmonds.pmi.no.reg
    score = bnlearn::score(net[[3]], categoric.dataset, type='loglik') 
    edmonds = list(pmi.no.reg.adj = adj.matrix.edmonds.pmi.no.reg, score = score, model=res)
    results[["edmonds"]] = edmonds

    # performs the reconstructions with Gabow
    res = NULL
    res = tronco.gabow(data,regularization=c("no_reg"), score=c('pmi'), epos = epos, eneg = eneg, silent = silent, boot.seed = seed)
    adj.matrix.gabow.pmi.no.reg = as.adj.matrix(res,model="gabow_no_reg_pmi")[['gabow_no_reg_pmi']]
    amat(net[[4]]) = adj.matrix.gabow.pmi.no.reg
    score = bnlearn::score(net[[4]], categoric.dataset, type='loglik') 
    gabow = list(pmi.no.reg.adj = adj.matrix.gabow.pmi.no.reg, score = score, model=res)
    results[["gabow"]] = gabow
     
    # performs the reconstructions with Chow Liu loglik, aic and bic
    res = NULL
    res = tronco.chowliu(data,regularization=c("loglik"), epos = epos, eneg = eneg, silent = silent, boot.seed = seed)   
    adj.matrix.chowliu.loglik = as.adj.matrix(res,model="chow_liu_loglik")[['chow_liu_loglik']]
    amat(net[[5]]) = adj.matrix.chowliu.loglik
    score = bnlearn::score(net[[5]], categoric.dataset, type='loglik')
    chowliu = list(loglik.adj=adj.matrix.chowliu.loglik, score = score, model=res)
    results[["chowliu"]] = chowliu
    
    # performs the reconstructions with Prim no_reg, loglik, aic and bic
    res = NULL
    res = tronco.prim(data,regularization=c("no_reg"), epos = epos, eneg = eneg, silent = silent, boot.seed = seed)
    adj.matrix.prim.no.reg = as.adj.matrix(res,model="prim_no_reg")[['prim_no_reg']]
    amat(net[[6]]) = adj.matrix.prim.no.reg
    score = bnlearn::score(net[[6]], categoric.dataset, type='loglik') 
    prim = list(no.reg.adj=adj.matrix.prim.no.reg, score = score, model=res)
    results[["prim"]] = prim
    
    return(results)
    
}

# perform the reconstructions
for(i in 1:nrow(dataset.missing.data)) {
    for (j in 1:ncol(dataset.missing.data)) {
        cat((((i - 1) * ncol(dataset.missing.data)) + j) , '/', nrow(dataset.missing.data) * ncol(dataset.missing.data), '\n')
        single.experiment = dataset.missing.data[[i,j]] 
        results = lapply(single.experiment, function(x, pass.error.rates){
            run.reconstructions(x$dataset,
                x$epos,
                x$eneg)
        })
        for (k in 1:length(results)) {
            dataset.missing.data[[i,j]][[k]]$reconstructions = results[[k]]
        }
    }
}

experiment.missing.data = dataset.missing.data
save(experiment.missing.data, file = "RData/experiment.missing.data.RData")
