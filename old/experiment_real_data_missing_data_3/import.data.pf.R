# load the required R packages
library(ape)
library(TRONCO)


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

colnames(original.data.paper3) = paste0('c', 1:ncol(original.data.paper3))
rownames(original.data.paper3) = paste0('s', 1:nrow(original.data.paper3))

stree = nj(dist.gene(original.data.paper3))
plot(stree)
dev.copy2pdf(file = 'data_3.pdf')
