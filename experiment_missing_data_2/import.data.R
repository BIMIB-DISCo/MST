# load the required R packages
library(mi)
library(TRONCO)

# set the seed
set.seed(12345)

# structure to save all the results
data_paper1 = list()

# read the data and format them
data_2 = read.table(file=paste0(getwd(),"/alterations.txt"),header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
rownames(data_2) = data_2[,1]
data_2 = data_2[,-1]
for (i in 2:nrow(data_2)) {
    for(j in 1:ncol(data_2)) {
        if(data_2[i,j]=="-") {
            data_2[i,j] = NA
        }
        else if(data_2[i,j]==data_2[1,j]) {
            data_2[i,j] = 0
        }
        else if(data_2[i,j]!=data_2[1,j]) {
            data_2[i,j] = 1
        }
    }
}
data_2 = data_2[-1,]
original.data.paper2 = data_2

# impute the missing data 100 times
dataset_imputed = mi(missing_data.frame(original.data.paper2))
dataset_imputed = complete(dataset_imputed,m=1)

# save the imputed datasets
data_imputations = list()
dataset.missing.data = matrix(list(), 1, 100)
colnames(dataset.missing.data) = paste("Experiment", 1:length(dataset_imputed))
rownames(dataset.missing.data) = "17"
for(i in 1:length(dataset_imputed)) {
    curr_imputation = dataset_imputed[[i]]
    curr_imputation = as.matrix(curr_imputation[,1:ncol(original.data.paper2)])
    curr_imputation = matrix(mapply(as.numeric,curr_imputation),nrow=nrow(original.data.paper2),ncol=ncol(original.data.paper2))
    colnames(curr_imputation) = colnames(original.data.paper2)
    rownames(curr_imputation) = rownames(original.data.paper2)
    exp = list()
    exp$dataset = curr_imputation
    exp$epos = 0
    exp$eneg = 0
    exp = list("1" = exp)
    dataset.missing.data[[1,i]] = exp
}

save(dataset.missing.data, file = "RData/dataset.missing.data.RData")
save(original.data.paper2, file = "RData/original.data.paper1.RData")


# create scite dataset 

scite = matrix(0, ncol = ncol(original.data.paper2), nrow = nrow(original.data.paper2))
colnames(scite) = colnames(original.data.paper2)
rownames(scite) = rownames(original.data.paper2)
for (i in 1:nrow(original.data.paper2)) {
    for (j in 1:ncol(original.data.paper2)) {
        if (is.na(original.data.paper2[[i,j]])) {
            scite[[i,j]] = 3
        } else if (original.data.paper2[[i,j]] == 1) {
            scite[[i,j]] = 1
        }
    }
}

save(scite, file='RData/scite.RData')

scite.dataset = NULL
scite.dataset$dataset = scite
scite.dataset$epos = 0
scite.dataset$eneg = 0

scite.input = matrix(list(), ncol = 1, nrow = 1)
colnames(scite.input) = 'Experiment 1'
rownames(scite.input) = '58'
scite.input[[1,1]] = list(scite.dataset)

source('../generate.scite.input.R')
create.scite.input(scite.input, 'single', 'missing', 0)
