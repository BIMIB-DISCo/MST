em.list = list()
em.list["2+3->5"] = 0
em.list["2-3->5"] = 0
em.list["3+4->6"] = 0
em.list["3-4->6"] = 0
em.list["2+3->5 + 3+4->6"] = 0
em.list["2-3->5 + 3-4->6"] = 0

count_caprese = em.list
count_capri = em.list
count_chow_liu = em.list
count_edmonds = em.list
count_gabow = em.list
count_prim = em.list
count_scite = em.list

count_conv = function(count.algo, adj.matrix) {
	if (adj.matrix[2,5] == 1
		&& adj.matrix[3,5] == 1
		&& adj.matrix[3,6] == 1
		&& adj.matrix[4,6] == 1) {
		print(adj.matrix)
		count.algo[["2+3->5 + 3+4->6"]] = count.algo[["2+3->5 + 3+4->6"]] + 1
	} else if ((adj.matrix[2,5] == 1
		|| adj.matrix[3,5] == 1)
		&& (adj.matrix[3,6] == 1
		|| adj.matrix[4,6] == 1)) {
		print(adj.matrix)
		count.algo[["2-3->5 + 3-4->6"]] = count.algo[["2-3->5 + 3-4->6"]] + 1
	} 
	if (adj.matrix[2,5] == 1
		&& adj.matrix[3,5] == 1){
		print(adj.matrix)
		count.algo[["2+3->5"]] = count.algo[["2+3->5"]] + 1
	} else if (adj.matrix[2,5] == 1
		|| adj.matrix[3,5] == 1){
		print(adj.matrix)
		count.algo[["2-3->5"]] = count.algo[["2-3->5"]] + 1
	} 
	if (adj.matrix[3,6] == 1
		&& adj.matrix[4,6] == 1) {
		print(adj.matrix)
		count.algo[["3+4->6"]] = count.algo[["3+4->6"]] + 1
	} else if (adj.matrix[3,6] == 1
		|| adj.matrix[4,6] == 1) {
		print(adj.matrix)
		count.algo[["3-4->6"]] = count.algo[["3-4->6"]] + 1
	}
	return(count.algo)
}

for (i in 1:100) {
	exp = experiments.multiple.biopses.convergent.scite[[1, i]][[1]]
	recon = exp$reconstruction
	
	caprese = recon$caprese$adj.caprese$caprese
	count_caprese = count_conv(count_caprese, caprese)
	
	capri = recon$capri$bic.adj$capri_bic
	count_capri = count_conv(count_capri, capri)

	chow_liu = recon$chowliu$loglik.adj$chow_liu_loglik
	count_chow_liu = count_conv(count_chow_liu, chow_liu)

	edmonds = recon$edmonds$pmi.no.reg.adj$edmonds_no_reg_pmi
	count_edmonds = count_conv(count_edmonds, edmonds)

	gabow = recon$gabow$pmi.no.reg.adj$gabow_no_reg_pmi
	count_gabow = count_conv(count_gabow, gabow)
	
	prim = recon$prim$no.reg.adj$prim_no_reg
	count_prim = count_conv(count_prim, prim)
	
	scite = recon$scite$no.reg.adj$scite_no_reg
	count_scite = count_conv(count_scite, scite)
	print(i)
}

