##############################################################################
###
### MST
###
### Import Data
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


## Load the required R packages

library(mi)
library(TRONCO)


## Set the seed

set.seed(333232)

epos = 1.24e-6
eneg = 0.0972


## Read the data and format them

data_1 = read.table(file = paste0(getwd(), "/alterations.txt"),
                    header = TRUE,
                    sep = "\t",
                    check.names = FALSE,
                    stringsAsFactors = FALSE)

for (i in 1:nrow(data_1)) {
    for(j in 1:ncol(data_1)) {
        if(data_1[i, j] == "-") {
            data_1[i, j] = NA
        }
        else if(data_1[i, j] == 2) {
            data_1[i, j] = 1
        }
    }
}

data_1 = as.matrix(data_1)
original.data.paper3 = apply(data_1, 2, as.numeric)
rownames(original.data.paper3) = 1:nrow(original.data.paper3)


## Impute the missing data 100 times

dataset_imputed = mi(missing_data.frame(original.data.paper3),
                     n.iter = 30,
                     n.chains = 4)
dataset_imputed = complete(dataset_imputed, m = 100)


## Save the imputed datasets

data_imputations = list()
dataset.missing.data = matrix(list(), 1, 100)
colnames(dataset.missing.data) = paste("Experiment", 1:length(dataset_imputed))
rownames(dataset.missing.data) = "47"
for(i in 1:length(dataset_imputed)) {
    curr_imputation = dataset_imputed[[i]]
    curr_imputation = as.matrix(curr_imputation[, 1:ncol(original.data.paper3)])
    curr_imputation = matrix(mapply(as.numeric, curr_imputation),
                             nrow = nrow(original.data.paper3),
                             ncol = ncol(original.data.paper3))
    colnames(curr_imputation) = colnames(original.data.paper3)
    rownames(curr_imputation) = rownames(original.data.paper3)
    curr_imputation[, colSums(curr_imputation, na.rm = TRUE) == 0] = 0
    exp = list()
    exp$dataset = curr_imputation
    exp$epos = epos
    exp$eneg = eneg
    exp = list("1" = exp)
    dataset.missing.data[[1, i]] = exp
}

save(dataset.missing.data, file = "RData/dataset.missing.data.RData")
save(original.data.paper3, file = "RData/original.data.paper3.RData")

### end of file -- import.data.R
