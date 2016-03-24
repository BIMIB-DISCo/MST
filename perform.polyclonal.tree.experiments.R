# simulate a dataset from single cells at given sample sizes and noise levels
"run.experiments.single.cells" = function ( type, true_tree, samples_num, nodes_probabilities = NA, e_pos, e_neg ) {
	
	# structure to save the results
	results = NULL
	
	# run the experiments
	if(type=="low") {
		for (i in samples_num) {
			for (j in 1:length(e_pos)) {
				curr_dataset = sample.single.cells.polyclonal.low(i,nodes_probabilities,e_pos[j],e_neg[j])
				results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
				results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
			}
		}
	}
	else if(type=="medium") {
		for (i in samples_num) {
			for (j in 1:length(e_pos)) {
				curr_dataset = sample.single.cells.polyclonal.medium(i,nodes_probabilities,e_pos[j],e_neg[j])
				results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
				results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
			}
		}
	}
	else if(type=="high") {
		for (i in samples_num) {
			for (j in 1:length(e_pos)) {
				curr_dataset = sample.single.cells.polyclonal.high(i,nodes_probabilities,e_pos[j],e_neg[j])
				results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
				results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
			}
		}
	}
	
	return(results)
	
}

# simulate a dataset from multiple biopses at given sample sizes and noise levels
"run.experiments.multiple.biopses" = function ( type, true_tree, samples_num, clones_probabilities, nodes_probabilities = NA, e_pos, e_neg, wild_type ) {
	
	
	# structure to save the results
	results = NULL
	
	# run the experiments
	if(type=="low") {
		for (i in samples_num) {
			for (j in 1:length(e_pos)) {
				curr_dataset = sample.multiple.biopses.polyclonal.low(i,clones_probabilities,nodes_probabilities,e_pos[j],e_neg[j],wild_type)
				results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
				results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
			}
		}
	}
	else if(type=="medium") {
		for (i in samples_num) {
			for (j in 1:length(e_pos)) {
				curr_dataset = sample.multiple.biopses.polyclonal.medium(i,clones_probabilities,nodes_probabilities,e_pos[j],e_neg[j],wild_type)
				results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
				results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
			}
		}
	}
	else if(type=="high") {
		for (i in samples_num) {
			for (j in 1:length(e_pos)) {
				curr_dataset = sample.multiple.biopses.polyclonal.high(i,clones_probabilities,nodes_probabilities,e_pos[j],e_neg[j],wild_type)
				results[[as.character(i)]][[as.character(j)]][["dataset"]] = curr_dataset
				results[[as.character(i)]][[as.character(j)]][["reconstructions"]] = run.reconstructions(curr_dataset,true_tree)
			}
		}
	}
	
	return(results)
	
}

# perform the reconstructions
"run.reconstructions" = function( dataset, true_tree ) {
	
	results = NULL
	
	# create a TRONCO object of the dataset
	data = import.genotypes(dataset)
	
	# performs the reconstructions with CAPRI loglik, aic and bic
	res = NULL
	res = tronco.capri(data,regularization=c("loglik","aic","bic"))
	adj.matrix.capri.loglik = as.adj.matrix(res,model="capri_loglik")
	results.capri.loglik = getStats(true_tree,adj.matrix.capri.loglik[["capri_loglik"]])
	adj.matrix.capri.aic = as.adj.matrix(res,model="capri_aic")
	results.capri.aic = getStats(true_tree,adj.matrix.capri.aic[["capri_aic"]])
	adj.matrix.capri.bic = as.adj.matrix(res,model="capri_bic")
	results.capri.bic = getStats(true_tree,adj.matrix.capri.bic[["capri_bic"]])
	capri = list(loglik.adj=adj.matrix.capri.loglik,loglik.res=results.capri.loglik,aic.adj=adj.matrix.capri.aic,aic.res=results.capri.aic,bic.adj=adj.matrix.capri.bic,bic.res=results.capri.bic)
	results[["capri"]] = capri
	
	# performs the reconstructions with CAPRESE
	res = NULL
	res = tronco.caprese(data)
	adj.matrix.caprese = as.adj.matrix(res,model="caprese")
	results.caprese = getStats(true_tree,adj.matrix.caprese[["caprese"]])
	caprese = list(adj.caprese=adj.matrix.caprese,caprese.res=results.caprese)
	results[["caprese"]] = caprese
	
	# performs the reconstructions with Edmonds no_reg, loglik, aic and bic
	res = NULL
	res = tronco.mst.edmonds(data,regularization=c("no_reg","loglik","aic","bic"))
	adj.matrix.edmonds.no.reg = as.adj.matrix(res,model="edmonds_no_reg")
	results.edmonds.no.reg = getStats(true_tree,adj.matrix.edmonds.no.reg[["edmonds_no_reg"]])	
	adj.matrix.edmonds.loglik = as.adj.matrix(res,model="edmonds_loglik")
	results.edmonds.loglik = getStats(true_tree,adj.matrix.edmonds.loglik[["edmonds_loglik"]])
	adj.matrix.edmonds.aic = as.adj.matrix(res,model="edmonds_aic")
	results.edmonds.aic = getStats(true_tree,adj.matrix.edmonds.aic[["edmonds_aic"]])
	adj.matrix.edmonds.bic = as.adj.matrix(res,model="edmonds_bic")
	results.edmonds.bic = getStats(true_tree,adj.matrix.edmonds.bic[["edmonds_bic"]])
	edmonds = list(no.reg.adj=adj.matrix.edmonds.no.reg,no.reg.res=results.edmonds.no.reg,loglik.adj=adj.matrix.edmonds.loglik,loglik.res=results.edmonds.loglik,aic.adj=adj.matrix.edmonds.aic,aic.res=results.edmonds.aic,bic.adj=adj.matrix.edmonds.bic,bic.res=results.edmonds.bic)
	results[["edmonds"]] = edmonds
	
	# performs the reconstructions with Chow Liu loglik, aic and bic
	res = NULL
	res = tronco.mst.chowliu(data,regularization=c("loglik","aic","bic"))	
	adj.matrix.chowliu.loglik = as.adj.matrix(res,model="chow_liu_loglik")
	results.chowliu.loglik = getStats(true_tree,adj.matrix.chowliu.loglik[["chow_liu_loglik"]])
	adj.matrix.chowliu.aic = as.adj.matrix(res,model="chow_liu_aic")
	results.chowliu.aic = getStats(true_tree,adj.matrix.chowliu.aic[["chow_liu_aic"]])
	adj.matrix.chowliu.bic = as.adj.matrix(res,model="chow_liu_bic")
	results.chowliu.bic = getStats(true_tree,adj.matrix.chowliu.bic[["chow_liu_bic"]])
	chowliu = list(loglik.adj=adj.matrix.chowliu.loglik,loglik.res=results.chowliu.loglik,aic.adj=adj.matrix.chowliu.aic,aic.res=results.chowliu.aic,bic.adj=adj.matrix.chowliu.bic,bic.res=results.chowliu.bic)
	results[["chowliu"]] = chowliu
	
	# performs the reconstructions with Prim no_reg, loglik, aic and bic
	res = NULL
	res = tronco.mst.prim(data,regularization=c("no_reg","loglik","aic","bic"))
	adj.matrix.prim.no.reg = as.adj.matrix(res,model="prim_no_reg")
	results.prim.no.reg = getStats(true_tree,adj.matrix.prim.no.reg[["prim_no_reg"]])
	adj.matrix.prim.loglik = as.adj.matrix(res,model="prim_loglik")
	results.prim.loglik = getStats(true_tree,adj.matrix.prim.loglik[["prim_loglik"]])
	adj.matrix.prim.aic = as.adj.matrix(res,model="prim_aic")
	results.prim.aic = getStats(true_tree,adj.matrix.prim.aic[["prim_aic"]])
	adj.matrix.prim.bic = as.adj.matrix(res,model="prim_bic")
	results.prim.bic = getStats(true_tree,adj.matrix.prim.bic[["prim_bic"]])
	prim = list(no.reg.adj=adj.matrix.prim.no.reg,no.reg.res=results.prim.no.reg,loglik.adj=adj.matrix.prim.loglik,loglik.res=results.prim.loglik,aic.adj=adj.matrix.prim.aic,aic.res=results.prim.aic,bic.adj=adj.matrix.prim.bic,bic.res=results.prim.bic)
	results[["prim"]] = prim
	
	return(results)
	
}
