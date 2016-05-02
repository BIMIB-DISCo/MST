library(TRONCO)
source('statistics.R')
source('edmonds.scores.R')

# set the true trees
low_true_tree = array(0,c(6,6))
low_true_tree[1,2] = 1
low_true_tree[1,3] = 1
low_true_tree[1,4] = 1
low_true_tree[2,5] = 1
low_true_tree[3,6] = 1

medium_true_tree = array(0,c(11,11))
medium_true_tree[1,2] = 1
medium_true_tree[1,3] = 1
medium_true_tree[1,4] = 1
medium_true_tree[2,5] = 1
medium_true_tree[4,6] = 1
medium_true_tree[5,7] = 1
medium_true_tree[5,8] = 1
medium_true_tree[6,9] = 1
medium_true_tree[6,10] = 1
medium_true_tree[6,11] = 1

high_true_tree = array(0,c(17,17))
high_true_tree[1,2] = 1
high_true_tree[1,3] = 1
high_true_tree[1,4] = 1
high_true_tree[2,5] = 1
high_true_tree[4,6] = 1
high_true_tree[5,7] = 1
high_true_tree[5,8] = 1
high_true_tree[6,9] = 1
high_true_tree[6,10] = 1
high_true_tree[6,11] = 1
high_true_tree[7,12] = 1
high_true_tree[9,13] = 1
high_true_tree[9,14] = 1
high_true_tree[9,15] = 1
high_true_tree[11,16] = 1
high_true_tree[11,17] = 1




compute.edmonds = function(dataset,
	sample.type,
	branching,
    true_tree) {

	e_pos_single_cells = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
	e_neg_single_cells = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
	sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
	e_pos_multiple_biopses = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
	e_neg_multiple_biopses = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
	sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)

	if (sample.type == "single") {
        epos = e_pos_single_cells
        eneg = e_neg_single_cells
        sample.size = sample_sizes_single_cells
    } else if (sample.type == "multiple") {
        epos = e_pos_multiple_biopses
        eneg = e_neg_multiple_biopses
        sample.size = sample_sizes_multiple_biopses
    } else {
        stop('sample.type must be "single" or "multiple"\n')
    }

    if (! branching %in% c('low', 'medium', 'high')) {
        stop('branching must be "low", "medium" or "high"')
    } 

    if (! dir.exists('RData')) {
    	stop('no buono')
    }
    count = 0
    for (sample in 1:nrow(dataset)) {
    	for (experiment in 1:ncol(dataset)) {
    		execution = dataset[[sample, experiment]]
    		for (noise in ls(execution)) {
    			noise_level = as.numeric(noise)
    			sample_size_level = sample.size[sample]
                object = execution[[noise]]
                data = object$dataset
                reconstructions = object$reconstructions
                # recon
                data = as.data.frame(data, StrongsAsFactors = FALSE)
                data = import.genotypes(data)
                # view(data)
                recon = tronco.mst.edmonds(data, silent = TRUE)

                score.standard = as.adj.matrix(recon)$edmonds_no_reg
                #print(score.standard)
                score.mi = perform.scores.edmonds(recon, 'MI')
                score.cmi = perform.scores.edmonds(recon, 'CMI')
                score.cmi2 = perform.scores.edmonds(recon, 'CMI2')
                score.cmi3 = perform.scores.edmonds(recon, 'CMI3')
                score.or = perform.scores.edmonds(recon, 'OR')
                score.or2 = perform.scores.edmonds(recon, 'OR2')
                #score.pearson = perform.scores.edmonds(recon, 'pearson')
                #score.kendall = perform.scores.edmonds(recon, 'kendall')
                #score.spearman = perform.scores.edmonds(recon, 'spearman')
                #score.causalor1 = perform.scores.edmonds(recon, 'CAUSAL_OR_1')
                #score.causalor2 = perform.scores.edmonds(recon, 'CAUSAL_OR_2')
                #score.causalor3 = perform.scores.edmonds(recon, 'CAUSAL_OR_3')
                #score.causalor4 = perform.scores.edmonds(recon, 'CAUSAL_OR_4')

                # adj matrix
                object$reconstructions$edmonds$no.reg.adj$edmonds_no_reg = score.standard
                object$reconstructions$edmonds$mi.adj$edmonds_no_reg = score.mi
                object$reconstructions$edmonds$cmi.adj$edmonds_no_reg = score.cmi
                object$reconstructions$edmonds$cmi2.adj$edmonds_no_reg = score.cmi2
                object$reconstructions$edmonds$cmi3.adj$edmonds_no_reg = score.cmi3
                object$reconstructions$edmonds$or.adj$edmonds_no_reg = score.or
                object$reconstructions$edmonds$or2.adj$edmonds_no_reg = score.or2
                #object$reconstructions$edmonds$pearson.adj$edmonds_no_reg = score.pearson
                #object$reconstructions$edmonds$kendall.adj$edmonds_no_reg = score.kendall
                #object$reconstructions$edmonds$spearman.adj$edmonds_no_reg = score.spearman
                #object$reconstructions$edmonds$causalor1.adj$edmonds_no_reg = score.causalor1
                #object$reconstructions$edmonds$causalor2.adj$edmonds_no_reg = score.causalor2
                #object$reconstructions$edmonds$causalor3.adj$edmonds_no_reg = score.causalor3
                #object$reconstructions$edmonds$causalor4.adj$edmonds_no_reg = score.causalor4

                statistics.standard = getStats(true_tree, score.standard)
                #print(statistics.standard)
                statistics.mi = getStats(true_tree, score.mi)
                statistics.cmi = getStats(true_tree, score.cmi)
                statistics.cmi2 = getStats(true_tree, score.cmi2)
                statistics.cmi3 = getStats(true_tree, score.cmi3)
                statistics.or = getStats(true_tree, score.or)
                statistics.or2 = getStats(true_tree, score.or2)
                #statistics.pearson = getStats(score.pearson, true_tree)
                #statistics.kendall = getStats(score.kendall, true_tree)
                #statistics.spearman = getStats(score.spearman, true_tree)
                #statistics.causalor1 = getStats(score.causalor1, true_tree)
                #statistics.causalor2 = getStats(score.causalor2, true_tree)
                #statistics.causalor3 = getStats(score.causalor3, true_tree)
                #statistics.causalor4 = getStats(score.causalor4, true_tree)

                object$reconstructions$edmonds$no.reg.res = statistics.standard
                object$reconstructions$edmonds$mi.res = statistics.mi
                object$reconstructions$edmonds$cmi.res = statistics.cmi
                object$reconstructions$edmonds$cmi2.res = statistics.cmi2
                object$reconstructions$edmonds$cmi3.res = statistics.cmi3
                object$reconstructions$edmonds$or.res = statistics.or
                object$reconstructions$edmonds$or2.res = statistics.or2
                #object$reconstructions$edmonds$pearson.res = statistics.pearson
                #object$reconstructions$edmonds$kendall.res = statistics.kendall
                #object$reconstructions$edmonds$spearman.res = statistics.spearman
                #object$reconstructions$edmonds$causalor1.res = statistics.causalor1
                #object$reconstructions$edmonds$causalor2.res = statistics.causalor2
                #object$reconstructions$edmonds$causalor3.res = statistics.causalor3
                #object$reconstructions$edmonds$causalor4.res = statistics.causalor4


                dataset[[sample, experiment]][[noise]] = object
                count = count+1
                print(count)
            }
        }
    }
    return(dataset)
}

