# load the required R packages
library(mi)
library(TRONCO)

# set the seed
set.seed(846232)

epos = 6.06e-6
eneg = 0.4309

# read the data and format them
data_1 = read.table(file=paste0(getwd(),"/alterations.txt"),header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
rownames(data_1) = data_1[,1]
data_1 = data_1[,-1]
for (i in 1:nrow(data_1)) {
    for(j in 1:ncol(data_1)) {
        if(data_1[i,j]=="-") {
            data_1[i,j] = NA
        }
        else if(data_1[i,j]==2) {
            data_1[i,j] = 1
        }
    }
}
origina.data.paper1 = data_1

# impute the missing data 100 times
dataset_imputed = mi(missing_data.frame(origina.data.paper1))
dataset_imputed = complete(dataset_imputed,m=100)

# save the imputed datasets
data_imputations = list()
dataset.missing.data = matrix(list(),1,100)
colnames(dataset.missing.data) = paste("Experiment",1:length(dataset_imputed))
rownames(dataset.missing.data) = "58"
for(i in 1:length(dataset_imputed)) {
    curr_imputation = dataset_imputed[[i]]
    curr_imputation = as.matrix(curr_imputation[,1:ncol(origina.data.paper1)])
    curr_imputation = matrix(mapply(as.numeric,curr_imputation),nrow=nrow(origina.data.paper1),ncol=ncol(origina.data.paper1))
    colnames(curr_imputation) = colnames(origina.data.paper1)
    rownames(curr_imputation) = rownames(origina.data.paper1)
    exp = list()
    exp$dataset = curr_imputation
    exp$epos = epos
    exp$eneg = eneg
    exp = list("1" = exp)
    dataset.missing.data[[1,i]] = exp
}

save(dataset.missing.data, file = "RData/dataset.missing.data.RData")
save(origina.data.paper1, file = "RData/original.data.paper1.RData")
