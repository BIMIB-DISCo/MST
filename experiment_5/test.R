## load library
#library(igraph)
#library(bnlearn)
#library(Rgraphviz)
#library(sna)
#library(oncoNEM)
#library(TRONCO)
#
## load the data
#load('RData_old/dataset.random.single.cells.5.nodes.RData')


reconstruct.and.plot <- function(this.experiment, scite.tree, exp.code) {
    epos = this.experiment$epos
    eneg = this.experiment$eneg
    true.tree = this.experiment$true_tree


    colnames(true.tree) = paste0('G', 1:ncol(true.tree))
    rownames(true.tree) = paste0('G', 1:ncol(true.tree))

    print(scite.tree)
    colnames(scite.tree) = paste0('G', 1:ncol(true.tree))
    rownames(scite.tree) = paste0('G', 1:ncol(true.tree))    

    dataset = this.experiment$dataset
    colnames(dataset) = paste0('G', 1:ncol(true.tree))
    categoric.dataset = data.frame(apply(dataset, 2, factor))
    for (name in colnames(categoric.dataset)) {
        levels(categoric.dataset[[name]]) = c(0,1)
    }

    #onconem.dataset = t(dataset)
    #
    #recon.onconem = oncoNEM$new(Data = onconem.dataset, FPR = epos, FNR = eneg)
    #recon.onconem$search(delta = 50, verbose = TRUE)
    ##linear.tree.onconem = recon.onconem$best$tree
    #score.onconem = recon.onconem$best$llh
    #oncoTree = clusterOncoNEM(oNEM = recon.onconem, epsilon = 10)
    #post = oncoNEMposteriors(tree = oncoTree$g,
    #    clones = oncoTree$clones,
    #    Data = recon.onconem$Data,
    #    FPR = recon.onconem$FPR,
    #    FNR = recon.onconem$FNR)
    #edgeLengths = colSums(post$p_theta)[-1]

    #dev.new()
    #plotTree(tree = oncoTree$g,
    #    clones = oncoTree$clones,
    #    e.length = edgeLengths,
    #    label.length = 4)
    #title(paste('\noncoNEM ', round(score.onconem, 3)))
    #dev.copy2pdf(file = paste0(exp.code, '_onconem.pdf'))
    #dev.off()


    #onconem.tree = matrix(0, ncol = length(linear.tree.onconem), nrow = length(linear.tree.onconem))
    #colnames(onconem.tree) = paste0('S', 1:length(linear.tree.onconem))
    #rownames(onconem.tree) = paste0('S', 1:length(linear.tree.onconem))
    #for (i in 1:length(linear.tree.onconem)) {
    #    onconem.tree[linear.tree.onconem[i], i] = 1
    #}




    #print(dataset)
    dataset = import.genotypes(dataset)
    print(dataset)
    #oncoprint(dataset)
    dev.new()
    par(mfrow =c(4, 4))

    true.tree.graph = graph.adjacency(true.tree)
    true.tree.nel = igraph.to.graphNEL(true.tree.graph)

    fontsize = rep(8, ncol(true.tree))
    names(fontsize) = paste0('G', 1:ncol(true.tree))

    plot(true.tree.nel, nodeAttrs = list(fontsize = fontsize))
    title('\nTRUE TREE')

    print(epos)
    print(eneg)

    recon.edmonds = tronco.mst.edmonds(dataset, 
        regularization = c("no_reg"),
        score = c('entropy', 'pmi', 'cpmi'),
        pvalue = 0.01,
        epos = epos,
        eneg = eneg)

    recon.gabow = tronco.mst.gabow(dataset, 
        regularization = c("no_reg"),
        score = c('entropy', 'pmi', 'cpmi', 'mi'),
        pvalue = 0.01,
        epos = epos,
        eneg = eneg)

    recon.prim = tronco.mst.prim(dataset, 
        regularization = c("no_reg"),
        pvalue = 0.01,
        epos = epos,
        eneg = eneg)

    recon.chowliu = tronco.mst.chowliu(dataset, 
        regularization = c("loglik"),
        pvalue = 0.01,
        epos = epos,
        eneg = eneg)

    recon.capri = tronco.capri(dataset, 
        regularization = c("bic"),
        pvalue = 0.01,
        epos = epos,
        eneg = eneg)

    recon.caprese = tronco.caprese(dataset, 
        epos = epos,
        eneg = eneg)

    recon.gabow.false = tronco.mst.gabow(dataset, 
        regularization = c("no_reg"),
        score = c('entropy', 'pmi', 'cpmi', 'mi'),
        pvalue = 0.01,
        epos = epos,
        eneg = eneg,
        do.raising = FALSE)

    #recon.mltree = tronco.mltree(dataset)

    #tree.list = recon.mltree$tree.list

    net = empty.graph(colnames(true.tree), num = 20)

    model.pmi.gabow = recon.gabow$model$gabow_no_reg_pmi$adj.matrix$adj.matrix.fit
    amat(net[[1]]) = model.pmi.gabow
    score.pmi.gabow = bnlearn::score(net[[1]], categoric.dataset, type='loglik')

    model.pmi.gabow.false = recon.gabow.false$model$gabow_no_reg_pmi$adj.matrix$adj.matrix.fit
    amat(net[[2]]) = model.pmi.gabow.false
    score.pmi.gabow.false = bnlearn::score(net[[2]], categoric.dataset, type='loglik')

    model.mi.gabow = recon.gabow$model$gabow_no_reg_mi$adj.matrix$adj.matrix.fit
    amat(net[[3]]) = model.mi.gabow
    score.mi.gabow = bnlearn::score(net[[3]], categoric.dataset, type='loglik')

    model.mi.gabow.false = recon.gabow.false$model$gabow_no_reg_mi$adj.matrix$adj.matrix.fit
    amat(net[[4]]) = model.mi.gabow.false
    score.mi.gabow.false = bnlearn::score(net[[4]], categoric.dataset, type='loglik')

    model.cpmi.gabow = recon.gabow$model$gabow_no_reg_cpmi$adj.matrix$adj.matrix.fit
    amat(net[[5]]) = model.cpmi.gabow
    score.cpmi.gabow = bnlearn::score(net[[5]], categoric.dataset, type='loglik')

    model.cpmi.gabow.false = recon.gabow.false$model$gabow_no_reg_cpmi$adj.matrix$adj.matrix.fit
    amat(net[[6]]) = model.cpmi.gabow.false
    score.cpmi.gabow.false = bnlearn::score(net[[6]], categoric.dataset, type='loglik')

    model.entropy.gabow = recon.gabow$model$gabow_no_reg_entropy$adj.matrix$adj.matrix.fit
    amat(net[[7]]) = model.entropy.gabow
    score.entropy.gabow = bnlearn::score(net[[7]], categoric.dataset, type='loglik')

    model.entropy.gabow.false = recon.gabow.false$model$gabow_no_reg_entropy$adj.matrix$adj.matrix.fit
    amat(net[[8]]) = model.entropy.gabow.false
    score.entropy.gabow.false = bnlearn::score(net[[8]], categoric.dataset, type='loglik')

    model.cpmi.edmonds = recon.edmonds$model$edmonds_no_reg_cpmi$adj.matrix$adj.matrix.fit
    amat(net[[9]]) = model.cpmi.edmonds
    score.edmonds.cpmi = bnlearn::score(net[[9]], categoric.dataset, type='loglik')

    model.pmi.edmonds = recon.edmonds$model$edmonds_no_reg_pmi$adj.matrix$adj.matrix.fit
    amat(net[[10]]) = model.pmi.edmonds
    score.edmonds.pmi = bnlearn::score(net[[10]], categoric.dataset, type='loglik')

    model.entropy.edmonds = recon.edmonds$model$edmonds_no_reg_entropy$adj.matrix$adj.matrix.fit
    amat(net[[11]]) = model.entropy.edmonds
    score.edmonds.entropy = bnlearn::score(net[[11]], categoric.dataset, type='loglik')

    model.bic.capri = recon.capri$model$capri_bic$adj.matrix$adj.matrix.fit
    amat(net[[12]]) = model.bic.capri
    score.capri.bic = bnlearn::score(net[[12]], categoric.dataset, type='loglik')

    model.caprese = recon.caprese$model$caprese$adj.matrix$adj.matrix.fit
    amat(net[[13]]) = model.caprese
    score.caprese = bnlearn::score(net[[13]], categoric.dataset, type='loglik')

    model.prim = recon.prim$model$prim_no_reg$adj.matrix$adj.matrix.fit
    amat(net[[14]]) = model.prim
    score.prim = bnlearn::score(net[[14]], categoric.dataset, type='loglik')

    model.chowliu = recon.chowliu$model$chow_liu_loglik$adj.matrix$adj.matrix.fit
    amat(net[[15]]) = model.chowliu
    score.chowliu = bnlearn::score(net[[15]], categoric.dataset, type='loglik')

    model.scite = scite.tree
    amat(net[[20]]) = model.scite
    score.scite = bnlearn::score(net[[20]], categoric.dataset, type='loglik')

    cat('scite \n')
    res.scite = getStats(true.tree, model.scite)

    #cat('pf cyclic \n')
    #print(pf.cyclic)
    #res.pf.cyclic = getStats(true.tree, pf.cyclic)

    #cat('pf \n')
    #res.pf = getStats(true.tree, pf)

    #pf.cyclic.graph = graph.adjacency(pf.cyclic)
    #pf.cyclic.nel = igraph.to.graphNEL(pf.cyclic.graph)
    #plot(pf.cyclic.nel)
    #title('\nPF CYCLIC')

    #pf.graph = graph.adjacency(pf)
    #pf.nel = igraph.to.graphNEL(pf.graph)
    #plot(pf.nel)
    #title('\nPF')

    #pmi.graph = graph.adjacency(model.pmi)
    #pmi.nel = igraph.to.graphNEL(pmi.graph)
    #plot(pmi.nel)
    #title(paste('\nPMI ', round(score.pmi, 3)))

    #cpmi.graph = graph.adjacency(model.cpmi)
    #cpmi.nel = igraph.to.graphNEL(cpmi.graph)
    #plot(cpmi.nel)
    #title(paste('\nCPMI', round(score.cpmi, 3)))

    #mi.graph.gabow = graph.adjacency(model.mi.gabow)
    #mi.nel.gabow = igraph.to.graphNEL(mi.graph.gabow)
    #plot(mi.nel.gabow, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nMI GABOW ', round(score.mi.gabow, 3)))

    pmi.graph.gabow = graph.adjacency(model.pmi.gabow)
    pmi.nel.gabow = igraph.to.graphNEL(pmi.graph.gabow)
    plot(pmi.nel.gabow, nodeAttrs = list(fontsize = fontsize))
    title(paste('\nPMI GABOW ', round(score.pmi.gabow, 3)))

    #cpmi.graph.gabow = graph.adjacency(model.cpmi.gabow)
    #cpmi.nel.gabow = igraph.to.graphNEL(cpmi.graph.gabow)
    #plot(cpmi.nel.gabow, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nCPMI GABOW ', round(score.cpmi.gabow, 3)))

    #entropy.graph.gabow = graph.adjacency(model.entropy.gabow)
    #entropy.nel.gabow = igraph.to.graphNEL(entropy.graph.gabow)
    #plot(entropy.nel.gabow, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nENTROPY GABOW ', round(score.entropy.gabow, 3)))


    #mi.graph.gabow.false = graph.adjacency(model.mi.gabow.false)
    #mi.nel.gabow.false = igraph.to.graphNEL(mi.graph.gabow.false)
    #plot(mi.nel.gabow.false, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nMI GABOW F ', round(score.mi.gabow.false, 3)))

    #pmi.graph.gabow.false = graph.adjacency(model.pmi.gabow.false)
    #pmi.nel.gabow.false = igraph.to.graphNEL(pmi.graph.gabow.false)
    #plot(pmi.nel.gabow.false, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nPMI GABOW F ', round(score.pmi.gabow.false, 3)))

    #cpmi.graph.gabow.false = graph.adjacency(model.cpmi.gabow.false)
    #cpmi.nel.gabow.false = igraph.to.graphNEL(cpmi.graph.gabow.false)
    #plot(cpmi.nel.gabow.false, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nCPMI GABOW F ', round(score.cpmi.gabow.false, 3)))

    #entropy.graph.gabow.false = graph.adjacency(model.entropy.gabow.false)
    #entropy.nel.gabow.false = igraph.to.graphNEL(entropy.graph.gabow.false)
    #plot(entropy.nel.gabow.false, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nENTROPY GABOW F', round(score.entropy.gabow.false, 3)))

    #cpmi.graph.edmonds = graph.adjacency(model.cpmi.edmonds)
    #cpmi.nel.edmonds = igraph.to.graphNEL(cpmi.graph.edmonds)
    #plot(cpmi.nel.edmonds, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nCPMI EDMONDS ', round(score.edmonds.cpmi, 3)))

    pmi.graph.edmonds = graph.adjacency(model.pmi.edmonds)
    pmi.nel.edmonds = igraph.to.graphNEL(pmi.graph.edmonds)
    plot(pmi.nel.edmonds, nodeAttrs = list(fontsize = fontsize))
    title(paste('\nPMI EDMONDS ', round(score.edmonds.pmi, 3)))

    #entropy.graph.edmonds = graph.adjacency(model.entropy.edmonds)
    #entropy.nel.edmonds = igraph.to.graphNEL(entropy.graph.edmonds)
    #plot(entropy.nel.edmonds, nodeAttrs = list(fontsize = fontsize))
    #title(paste('\nENTROPY EDMONDS ', round(score.edmonds.entropy, 3)))

    bic.graph.edmonds = graph.adjacency(model.bic.capri)
    bic.nel.edmonds = igraph.to.graphNEL(bic.graph.edmonds)
    plot(bic.nel.edmonds, nodeAttrs = list(fontsize = fontsize))
    title(paste('\nCAPRI BIC ', round(score.capri.bic, 3)))

    graph.caprese = graph.adjacency(model.caprese)
    nel.caprese = igraph.to.graphNEL(graph.caprese)
    plot(nel.caprese, nodeAttrs = list(fontsize = fontsize))
    title(paste('\nCAPRESE ', round(score.caprese, 3)))

    graph.prim = graph.adjacency(model.prim)
    nel.prim = igraph.to.graphNEL(graph.prim)
    plot(nel.prim, nodeAttrs = list(fontsize = fontsize))
    title(paste('\nPRIM ', round(score.prim, 3)))

    graph.chowliu = graph.adjacency(model.chowliu)
    nel.chowliu = igraph.to.graphNEL(graph.chowliu)
    plot(nel.chowliu, nodeAttrs = list(fontsize = fontsize))
    title(paste('\nCHOW LIU ', round(score.chowliu, 3)))

    #title(paste('\nMLTREE ', round(score.mltree, 3)))

    scite.graph = graph.adjacency(model.scite)
    scite.nel = igraph.to.graphNEL(scite.graph)
    plot(scite.nel, nodeAttrs = list(fontsize = fontsize))
    title(paste('\nSCITE ', round(score.scite, 3)))

    if (FALSE && ((res.pf$tp != res.pf.cyclic$tp)
        || (res.pf$tn != res.pf.cyclic$tn))) {

        cat('TRUE TREE\n')
        print(true.tree)
        cat('PF CYCLIC\n')
        print(pf.cyclic)
        cat('PF\n')
        print(pf)
        cat('SCITE\n')
        print(model.scite)

        cat ("Press [enter] to continue")
        line = readline()
    }

    #cat('Type "S" to save output: ')
    #line = readline()
    #if (line %in% c('s', 'S')) {

    if (TRUE || !all(true.tree == model.scite)
        || !all(true.tree == model.pmi.gabow.false)) {
        #|| res.pf$fn > 0) {

        dev.copy2pdf(file = paste0(exp.code, '_score.pdf'))
        dev.off()

        print(exp.code)
        dev.new()
        oncoprint(recon.gabow, title = exp.code, file = paste0(exp.code, '_oncoprint.pdf'))
        dev.off()

        #dev.new()
        #par(mfrow =c(5, 5))
        #for (tree in tree.list) {
        #    tree.graph = graph.adjacency(tree)
        #    tree.nel = igraph.to.graphNEL(tree.graph)
        #    net = empty.graph(colnames(true.tree))
        #    amat(net) = tree
        #    score = bnlearn::score(net, categoric.dataset, type='loglik')
        #    #cat('score: ', score, '\n')
        #    #bnlearn:::print.bn(net)
        #    plot(tree.nel)
        #    title(paste('\nMLTREE ', round(score, 3)))
        #}
        #dev.copy2pdf(file = paste0(exp.code, '_mltree.pdf'))
        #dev.off()

    } else {
        dev.off()
    }

    #cat('marginal probs\n')
    #print(recon$model$mltree_no_reg$probabilities$probabilities.observed$marginal.probs)
    #cat('joint probs\n')
    #print(recon$model$mltree_no_reg$probabilities$probabilities.observed$joint.probs)
    #cat('conditional probs\n')
    #print(recon$model$mltree_no_reg$probabilities$probabilities.observed$conditional.probs)
    #cat('selective advantage relations pf\n')
    ##print(as.selective.advantage.relations(recon.mltree, type='pf')[['pf']])

    return(recon.gabow.false)
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
        for (j in 1:ncol(inferred_matrix)) {
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

    cat('tp: ', tp, ' tn: ', tn, ' fp: ', fp, ' fn: ', fn, '\n')

    # return the results
    results_values = list(accuracy=accuracy,sensitivity=sensitivity,specificity=specificity,hamming_distance=hamming_distance,
        tp = tp, tn = tn, fp = fp, fn = fn)
    return(results_values)
    
}



# 3 1 38


# select a dataset
sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
e_pos_single_cells = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
e_neg_single_cells = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
noise = 3
sample = 5
sample.scite = sample_sizes_single_cells[[sample]]
epos = e_pos_single_cells[[noise]]
eneg = e_neg_single_cells[[noise]]
branching = 'high'
sample.type = 'single'
#dev.new()

dataset = dataset.single.cells.random.columns.high

#for (exp in 1:ncol(dataset.random.single.cells.5.nodes)) {
for (exp in 1:100) {
    print(exp)
    exp.code = paste(branching, sample.type, sample, exp, noise, sep = '_')
    exp.code = paste0('test/', exp.code)
    this.experiment = dataset[[sample, exp]][[noise]]

    filename = paste0('scite_output/datasets/', branching,  '/', sample.type, '/', sample.scite, '_', exp, '_', noise, '_ml0')
    readdot = read.dot(paste0(filename, '.correct.gv'))
    graph = graph.adjacency(readdot)
    scite.tree.raw = get.adjacency(graph, sparse = FALSE)

    scite.tree = matrix(0, ncol = ncol(scite.tree.raw) - 1, nrow = ncol(scite.tree.raw) - 1)
    for (i in 1:nrow(scite.tree)) {
        for(j in 1:ncol(scite.tree)) {
            if (scite.tree.raw[as.character(i), as.character(j)] == 1) {
                scite.tree[i,j] = 1
            }
        }
    }

    reconstruct.and.plot(this.experiment, scite.tree, exp.code)
}






#dev.new()
#par(mfrow =c(1, 2))
#this.experiment = dataset.random.single.cells.5.nodes[[1, 38]][[3]]
#data = this.experiment$dataset
#true.tree = this.experiment$true_tree$structure
#net = empty.graph(colnames(true.tree), num = 2)
#amat(net[[1]]) = true.tree
#
#test.graph = true.tree
#test.graph[[1,3]] = 0
#test.graph[[1,4]] = 1
#test.graph[[4,3]] = 1
#test.graph[[3,4]] = 0
#
#amat(net[[2]]) = test.graph
#scores = rep(0, 2)
#data = data.frame(apply(data, 2, factor))
#scores[[1]] = score(net[[1]], data, type='loglik')
#scores[[2]] = score(net[[2]], data, type='loglik')
#scores
#
#
#tt.graph = graph.adjacency(true.tree)
#tt.nel = igraph.to.graphNEL(tt.graph)
#plot(tt.nel)
#title(paste('\nTT', scores[[1]]))
#
#err.graph = graph.adjacency(test.graph)
#err.nel = igraph.to.graphNEL(err.graph)
#plot(err.nel)
#title(paste('\nERR', scores[[2]]))