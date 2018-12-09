##############################################################################
###
### MST
###
### Generate Main
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


load('RData/dataset.single.cells.medium.RData')
source('../generate.scite.input.R')


dataset = dataset.single.cells.medium['50', , drop = F]
mini_dataset_01 = dataset

epos1 = c(0.003, 0.004, 0.005, 0.006, 0.007)
eneg1 = rep(0.03, 5)

for (i in 1:ncol(dataset)) {
    this_dataset = dataset[[1, i]][[2]]
    new_list_of_dataset = list()
    for (j in 1:length(epos1)) {
        this_dataset$epos = epos1[[j]]
        this_dataset$eneg = eneg1[[j]]
        new_list_of_dataset[[as.character(j)]] = this_dataset
    }
    mini_dataset_01[[1, i]] = new_list_of_dataset
}

save(mini_dataset_01, file = 'RData/mini_dataset_01.RData')
create.scite.input(mini_dataset_01, 'single', 'mini_01', scite.sd)


dataset = dataset.single.cells.medium['50', , drop = F]
mini_dataset_02 = dataset

epos2 = c(0.003, 0.004, 0.005, 0.006, 0.007)
eneg2 = rep(0.04, 5)

for (i in 1:ncol(dataset)) {
    this_dataset = dataset[[1, i]][[2]]
    new_list_of_dataset = list()
    for (j in 1:length(epos2)) {
        this_dataset$epos = epos2[[j]]
        this_dataset$eneg = eneg2[[j]]
        new_list_of_dataset[[as.character(j)]] = this_dataset
    }
    mini_dataset_02[[1, i]] = new_list_of_dataset
}

save(mini_dataset_02, file = 'RData/mini_dataset_02.RData')
create.scite.input(mini_dataset_02, 'single', 'mini_02', scite.sd)


dataset = dataset.single.cells.medium['50', , drop = F]
mini_dataset_03 = dataset

epos3 = c(0.003, 0.004, 0.005, 0.006, 0.007)
eneg3 = rep(0.05, 5)

for (i in 1:ncol(dataset)) {
    this_dataset = dataset[[1, i]][[2]]
    new_list_of_dataset = list()
    for (j in 1:length(epos3)) {
        this_dataset$epos = epos3[[j]]
        this_dataset$eneg = eneg3[[j]]
        new_list_of_dataset[[as.character(j)]] = this_dataset
    }
    mini_dataset_03[[1, i]] = new_list_of_dataset
}

save(mini_dataset_03, file = 'RData/mini_dataset_03.RData')
create.scite.input(mini_dataset_03, 'single', 'mini_03', scite.sd)


dataset = dataset.single.cells.medium['50', , drop = F]
mini_dataset_04 = dataset

epos4 = c(0.003, 0.004, 0.005, 0.006, 0.007)
eneg4 = rep(0.06, 5)

for (i in 1:ncol(dataset)) {
    this_dataset = dataset[[1, i]][[2]]
    new_list_of_dataset = list()
    for (j in 1:length(epos4)) {
        this_dataset$epos = epos4[[j]]
        this_dataset$eneg = eneg4[[j]]
        new_list_of_dataset[[as.character(j)]] = this_dataset
    }
    mini_dataset_04[[1, i]] = new_list_of_dataset
}

save(mini_dataset_04, file = 'RData/mini_dataset_04.RData')
create.scite.input(mini_dataset_04, 'single', 'mini_04', scite.sd)

dataset = dataset.single.cells.medium['50', , drop = F]
mini_dataset_05 = dataset

epos5 = c(0.003, 0.004, 0.005, 0.006, 0.007)
eneg5 = rep(0.07, 5)

for (i in 1:ncol(dataset)) {
    this_dataset = dataset[[1, i]][[2]]
    new_list_of_dataset = list()
    for (j in 1:length(epos5)) {
        this_dataset$epos = epos5[[j]]
        this_dataset$eneg = eneg5[[j]]
        new_list_of_dataset[[as.character(j)]] = this_dataset
    }
    mini_dataset_05[[1, i]] = new_list_of_dataset
}

save(mini_dataset_05, file = 'RData/mini_dataset_05.RData')
create.scite.input(mini_dataset_05, 'single', 'mini_05', scite.sd)

### end of file -- generate.main.R
