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

library(mi)


## Set the seed

set.seed(12345)

source('../generate.scite.input.R')
missing.levels = c(10, 20, 30, 40)


load('RData/dataset.single.cells.medium.RData')
dataset.missing = dataset.single.cells.medium['75', 1:10, drop = F]

for (i in 1:nrow(dataset.missing)) {
    for (j in 1:ncol(dataset.missing)) {
        data = dataset.missing[[i, j]]
        noise.free = data[[1]]
        missing.data.dataset = list()
        for (missing.level in missing.levels) {
            dataset = noise.free
            for (k in 1:as.integer(length(dataset$dataset) / 100 * missing.level)) {
                sample = sample(1:nrow(dataset$dataset), 1)
                gene = sample(1:ncol(dataset$dataset), 1)
                while (is.na(dataset$dataset[[sample, gene]])) {
                    sample = sample(1:nrow(dataset$dataset), 1)
                    gene = sample(1:ncol(dataset$dataset), 1)
                }
                dataset$dataset[[sample, gene]] = NA
            }
            missing.data.dataset[[as.character(k)]] = dataset
        }
        dataset.missing[[i, j]] = missing.data.dataset
    }
}

cat('scite\n')
create.scite.input(dataset.missing, 'single', 'missing', scite.sd)

dataset.missing.impute = matrix(list(), nrow = 4, ncol = 10)

for (missing.level in 1:nrow(dataset.missing.impute)) {
    for (exp in 1:ncol(dataset.missing.impute)) {
        data = dataset.missing[[1, exp]][[missing.level]]
        dataset.imputed = mi(missing_data.frame(data$dataset))
        dataset.imputed = complete(dataset.imputed, m = 100)
        imputated.datasets = list()
        for(i in 1:length(dataset.imputed)) {
            curr.imputation = dataset.imputed[[i]]
            curr.imputation = as.matrix(curr.imputation[, 1:ncol(data$dataset)])
            curr.imputation = matrix(mapply(as.numeric, curr.imputation),
                                     nrow = nrow(data$dataset),
                                     ncol = ncol(data$dataset))
            colnames(curr.imputation) = colnames(data$dataset)
            rownames(curr.imputation) = rownames(data$dataset)
            curr.imputation[, colSums(curr.imputation, na.rm = TRUE) == 0] = 0
            imp.dataset = list()
            imp.dataset$dataset = curr.imputation
            imp.dataset$epos = data$epos
            imp.dataset$eneg = data$eneg
            imp.dataset$true_tree = data$true_tree
            imputated.datasets[[i]] = imp.dataset
        }
        dataset.missing.impute[[missing.level, exp]] = imputated.datasets
    }
}

rownames(dataset.missing.impute) = missing.levels
colnames(dataset.missing.impute) = paste('Experiment', 1:10)
save(dataset.missing.impute, file = 'RData/dataset.missing.impute.RData')

### end of file -- generate.main.R
