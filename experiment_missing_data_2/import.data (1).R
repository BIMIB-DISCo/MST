

# load the required R packages
library(mi)
library(TRONCO)

# set the seed
set.seed(12345)

# structure to save all the results
data_paper2 = list()

# read the data and format them
data_2 = read.table(file=paste0("alterations.txt"),header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
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
data_paper2[["original"]] = data_2


# library(pheatmap)
# t = data_paper2[["original"]]
# t[is.na(t)] = -1
# pheatmap(as.matrix(t))

# impute the missing data 100 times
dataset_imputed = mi(missing_data.frame(data_paper2[["original"]]))
dataset_imputed = complete(dataset_imputed,m=1)

# save the imputed datasets
data_imputations = list()
data_tronco = list()
for(i in 1:length(dataset_imputed)) {
    curr_imputation = dataset_imputed[[i]]
    curr_imputation = as.matrix(curr_imputation[,1:ncol(data_paper2[["original"]])])
    curr_imputation = matrix(mapply(as.numeric,curr_imputation),nrow=nrow(data_paper2[["original"]]),ncol=ncol(data_paper2[["original"]]))
    colnames(curr_imputation) = colnames(data_paper2[["original"]])
    rownames(curr_imputation) = rownames(data_paper2[["original"]])
    data_imputations[[i]] = curr_imputation
    data_tronco[[i]] = import.genotypes(curr_imputation)
}

# save the results
data_paper2[["imputations"]] = data_imputations
data_paper2[["tronco"]] = data_tronco
save(data_paper2,file="data_paper2.RData")
