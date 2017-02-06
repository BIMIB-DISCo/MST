# load the required R packages
library(TRONCO)
library(ape)

# set the seed

# structure to save all the results
data_paper1 = list()

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

colnames(origina.data.paper1) = paste0('c', 1:ncol(origina.data.paper1))
rownames(origina.data.paper1) = paste0('s', 1:nrow(origina.data.paper1))

stree = nj(dist.gene(origina.data.paper1))
plot(stree)
dev.copy2pdf(file = 'data_1.pdf')
