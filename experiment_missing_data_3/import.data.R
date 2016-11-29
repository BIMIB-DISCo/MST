# load the required R packages
library(mi)
library(TRONCO)

# set the seed
set.seed(343434)

epos = 1.24e-6
eneg = 0.0972

# structure to save all the results
data_paper1 = list()

# read the data and format them
data_1 = read.table(file=paste0(getwd(),"/alterations.txt"),header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
#rownames(data_1) = data_1[,1]
#data_1 = data_1[,-1]
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
data_1 = as.matrix(data_1)
original.data.paper3 = apply(data_1, 2, as.numeric)
rownames(original.data.paper3) = 1:nrow(original.data.paper3)



# impute the missing data 100 times
dataset_imputed = mi(missing_data.frame(original.data.paper3), n.iter = 30, n.chains = 4)
dataset_imputed = complete(dataset_imputed,m=100)

# save the imputed datasets
data_imputations = list()
dataset.missing.data = matrix(list(), 1, 100)
colnames(dataset.missing.data) = paste("Experiment", 1:length(dataset_imputed))
rownames(dataset.missing.data) = "47"
for(i in 1:length(dataset_imputed)) {
    curr_imputation = dataset_imputed[[i]]
    curr_imputation = as.matrix(curr_imputation[,1:ncol(original.data.paper3)])
    curr_imputation = matrix(mapply(as.numeric,curr_imputation),nrow=nrow(original.data.paper3),ncol=ncol(original.data.paper3))
    colnames(curr_imputation) = colnames(original.data.paper3)
    rownames(curr_imputation) = rownames(original.data.paper3)
    curr_imputation[,colSums(curr_imputation, na.rm = TRUE) == 0] = 0
    exp = list()
    exp$dataset = curr_imputation
    exp$epos = epos
    exp$eneg = eneg
    exp = list("1" = exp)
    dataset.missing.data[[1,i]] = exp
}

save(dataset.missing.data, file = "RData/dataset.missing.data.RData")
save(original.data.paper3, file = "RData/original.data.paper3.RData")


# create scite dataset 

scite = matrix(0, ncol = ncol(original.data.paper3), nrow = nrow(original.data.paper3))
colnames(scite) = colnames(original.data.paper3)
rownames(scite) = rownames(original.data.paper3)
for (i in 1:nrow(original.data.paper3)) {
    for (j in 1:ncol(original.data.paper3)) {
        if (is.na(original.data.paper3[[i,j]])) {
            scite[[i,j]] = 3
        } else if (original.data.paper3[[i,j]] == 1) {
            scite[[i,j]] = 1
        }
    }
}

save(scite, file='RData/scite.RData')

scite.dataset = NULL
scite.dataset$dataset = scite
scite.dataset$epos = epos
scite.dataset$eneg = eneg

scite.input = matrix(list(), ncol = 1, nrow = 1)
colnames(scite.input) = 'Experiment 1'
rownames(scite.input) = '47'
scite.input[[1,1]] = list(scite.dataset)

source('../generate.scite.input.R')
create.scite.input(scite.input, 'single', 'missing', 0)



new.names = paste0(as.character(round(colSums(scite == 1) / (rep(nrow(scite), ncol(scite)) - colSums(scite == 3)),3)*100), '%')
scite.new = scite
colnames(scite.new) = paste(colnames(scite), new.names)
pheatmap(t(scite.new), cluster_rows=F, cluster_cols=F, show_colnames=F, color=c('lightgrey','lightgreen','white', 'red'), legend = F, border_color='white')
dev.copy2pdf(file = 'plot/ph.pdf')