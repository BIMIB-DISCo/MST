# data - dummy
data = matrix(sample(0:1, 333*6, replace=TRUE),333,6)
data = data.frame(apply(data, 2, factor))
colnames(data) = as.character(1:6)
head(data)

# prima facie -- example
pf = matrix(0, 6, 6)
pf[1,1:4] = 1
pf[3,2:4] = 1
pf[3,6] = 1
pf[4,2] = 1
pf[5,1] = 1
diag(pf) = 0
pf

library(igraph)
G = graph.adjacency(pf)
plot(G)


########### Algo
# code generation
mywhich = function(x) {
    if(sum(pf[,x]) > 0) return(paste('which(pf[,', x,'] == 1)'))
    else return('0')
}

# Generate a code to exploit this command
# expand.grid(which(pf[,1] == 1), ...., which(pf[,nrow(pf)] == 1))
com = paste(sapply(1:nrow(pf), mywhich))
com = paste(com, collapse= ',')
com = paste('expand.grid(', com,')', collapse = '')
indexes = eval(parse(text=com))
print(paste(nrow(indexes), 'trees will be created.'))

# Generate all matrices - amazing
candidates = Reduce(
    append,
    apply(indexes, 1, function(x) {
        M = matrix(0L, nrow=ncol(pf), ncol = ncol(pf))
        for(i in 1:nrow(M)) M[x[i], i] = 1L
        return(list(M))
    })
)
print(paste("Total number of prima facie trees: ", length(candidates), sep = ''))

candidates
pf

# We go to bnlearn
library(bnlearn)
MAXRESULTS = 10 # output results -- full posterior by setting it to NUMTESTS
NUMTESTS = length(candidates)

# empty networks
net = empty.graph(as.character(1:ncol(pf)), num = NUMTESTS)

# set edges/compute score
scores = rep(0, NUMTESTS)
for(i in 1:NUMTESTS) # make this parallel
{
        amat(net[[i]]) = candidates[[i]]
        scores[i] = score(net[[i]], data, type = "loglik")
}

# sort and trim results
net = net[order(scores)]
net = net[1: min(MAXRESULTS, length(net))]
scores = sort(scores)

## extras for visualization
sq = ceiling(sqrt(length(net)))

plot(G)
title("ADJ")

dev.new()
par(mfrow =c(sq, sq))
for(i in 1: length(net))
{
        graphviz.plot(net[[i]])
        title(scores[i])
}