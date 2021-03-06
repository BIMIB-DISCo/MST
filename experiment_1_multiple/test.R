##############################################################################
###
### MST
###
### Test
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################

em.list = list()
em.list[1] = 0
em.list[2] = 0
em.list[3] = 0
em.list[4] = 0
em.list[5] = 0
em.list[6] = 0
count_caprese = em.list
count_capri = em.list
count_chow_liu = em.list
count_edmonds = em.list
count_gabow = em.list
count_prim = em.list
count_scite = em.list

for (i in 1:100) {
    exp = experiments.multiple.biopses.random.forest.scite[[1, i]][[1]]
    recon = exp$reconstruction
    
    ## caprese
    caprese = recon$caprese$adj.caprese$caprese
    gr = NULL
    gr = graph.adjacency(caprese)
    count_caprese[[count_components(gr)]] =
        count_caprese[[count_components(gr)]] + 1

    ## capri
    capri = recon$capri$bic.adj$capri_bic
    gr = NULL
    gr = graph.adjacency(capri)
    count_capri[[count_components(gr)]] =
        count_capri[[count_components(gr)]] + 1

    ## chow_liu
    chow_liu = recon$chowliu$loglik.adj$chow_liu_loglik
    gr = NULL
    gr = graph.adjacency(chow_liu)
    count_chow_liu[[count_components(gr)]] =
        count_chow_liu[[count_components(gr)]] + 1

    ## edmonds
    edmonds = recon$edmonds$pmi.no.reg.adj$edmonds_no_reg_pmi
    gr = NULL
    gr = graph.adjacency(edmonds)
    count_edmonds[[count_components(gr)]] =
        count_edmonds[[count_components(gr)]] + 1

    ## gabow
    gabow = recon$gabow$pmi.no.reg.adj$gabow_no_reg_pmi
    gr = NULL
    gr = graph.adjacency(gabow)
    count_gabow[[count_components(gr)]] =
        count_gabow[[count_components(gr)]] + 1

    ## prim
    prim = recon$prim$no.reg.adj$prim_no_reg
    gr = NULL
    gr = graph.adjacency(prim)
    count_prim[[count_components(gr)]] =
        count_prim[[count_components(gr)]] + 1

    ## scite
    scite = recon$scite$no.reg.adj$scite_no_reg
    gr = NULL
    gr = graph.adjacency(scite)
    count_scite[[count_components(gr)]] =
        count_scite[[count_components(gr)]] + 1
}

### end of file -- test.R
