

for (sample in rownames(dataset)) {
	for (exp in colnames(dataset)) {
		for (noise in names(dataset[[sample, exp]])) {
			cfg = dataset[[sample, exp]][[noise]]

			pmi.adj = cfg$reconstructions$edmonds$pmi.no.reg.adj$edmonds_no_reg_pmi
			cpmi.adj = cfg$reconstructions$edmonds$cpmi.no.reg.adj$edmonds_no_reg_cpmi
			true.tree = cfg$true_tree$structure

			results.edmonds.pmi.no.reg = getStats(true.tree, pmi.adj)
			results.edmonds.cpmi.no.reg = getStats(true.tree, cpmi.adj)

			cfg$reconstructions$edmonds$pmi.no.reg.res = results.edmonds.pmi.no.reg
			cfg$reconstructions$edmonds$cpmi.no.reg.res = results.edmonds.cpmi.no.reg


			dataset[[sample, exp]][[noise]] = cfg
		}
	}
}
