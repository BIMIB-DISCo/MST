#load('RData_old/dataset.random.single.cells.5.nodes.RData')

e_pos_single_cells = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035)
e_neg_single_cells = c(0.000, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350)
noise = 3
sample = 5
epos = e_pos_single_cells[[noise]]
eneg = e_neg_single_cells[[noise]]
exp = 13
nboot = 1000

this.experiment = dataset.random.single.cells.5.nodes[[sample, exp]][[noise]]

dataset = this.experiment$dataset

dataset = import.genotypes(dataset)

#recon = tronco.mst.edmonds(dataset, 
#    regularization = c("no_reg"),
#    score = c('entropy', 'pmi', 'cpmi'),
#    nboot = nboot,
#    epos = epos,
#    eneg = eneg)


recon = tronco.mst.gabow(dataset, 
    regularization = c("no_reg"),
    score = c('entropy', 'pmi', 'cpmi'),
    nboot = nboot,
    epos = epos,
    eneg = eneg,
    do.raising = TRUE)


load('scores.RData')

#dev.new()
#par(mfrow =c(1, 2))

from = 2
to = 4


# from -> to  
# a -> b
# 2 -> 4


# p b 
p.from = scores$marginal.probs.distributions[[from,1]]
#p.from = mean(p.from)

# p a 
p.to = scores$marginal.probs.distributions[[to,1]]
#p.to = mean(p.to)

# p a,b 
p.joint.from.to = scores$joint.probs.distributions[[from, to]]


# p(b|a)  = p(a,b) / p(a) 
p.to.given.from = p.joint.from.to / p.from

# p(b|-a)  =  p(b) - p(a,b) / 1 - p(a)
p.to.not.given.from = (p.to - p.joint.from.to) / (1 - p.from)

tp = recon$confidence[1,][[1]][[from, to]]
pr = recon$confidence[2,][[1]][[from, to]]
hg = recon$confidence[3,][[1]][[from, to]]

density.value = cbind(p.from, p.to)
colnames(density.value) = c(paste0('p(', from, ')'), paste0('p(', to, ')'))

density.value = as.data.frame(density.value)
density.value.s = stack(density.value)

density.raising.value = cbind(p.to.given.from, p.to.not.given.from)
colnames(density.raising.value) = c(paste0('p(', to, '|', from, ')'), paste0('p(', to, '|-', from, ')'))

density.raising.value = as.data.frame(density.raising.value)
density.raising.value.s = stack(density.raising.value)

hmax = 4

p1 = ggplot(density.value.s, aes(x = values)) + 
    geom_histogram(aes(y = ..density.., group = ind, colour=ind), position="identity", size = 0, fill = NA) + 
    geom_density(aes(group = ind, colour=ind), alpha=.5) +
    ggtitle("Temporal Priority") + 
    geom_text(aes(label = paste('nboot =', nboot), x = Inf, y = hmax), hjust = 0) +
    geom_text(aes(label = paste('epos =', epos), x = Inf, y = hmax - 0.5), hjust = 0) +
    geom_text(aes(label = paste('eneg =', eneg), x = Inf, y = hmax - 1), hjust = 0) +
    geom_text(aes(label = paste('tp =', tp), x = Inf, y = hmax - 1.5), hjust = 0) +
    geom_text(aes(label = paste('pr =', pr), x = Inf, y = hmax - 2), hjust = 0) +
    geom_text(aes(label = paste('hg =', hg), x = Inf, y = hmax - 2.5), hjust = 0) +
    geom_text(aes(label = paste0('p(',from , ') =', mean(p.from)), x = Inf, y = hmax - 3), hjust = 0) +
    geom_text(aes(label = paste0('p(',to , ') =', mean(p.to)), x = Inf, y = hmax - 3.5), hjust = 0) + 
    geom_text(aes(label = paste0('p(',from, ',', to , ') =', mean(p.joint.from.to)), x = Inf, y = hmax - 4), hjust = 0)
    
p2 = ggplot(density.raising.value.s, aes(x=values)) + 
    geom_histogram(aes(y = ..density.., group = ind, colour=ind), position="identity", size = 0, fill = NA) + 
    geom_density(aes(group=ind, colour=ind), alpha=.5) +
    ggtitle("Probability Raising") 



p1 = ggplot_gtable(ggplot_build(p1))
p1$layout$clip[p1$layout$name=="panel"] <- "off"

p2 = ggplot_gtable(ggplot_build(p2))
p2$layout$clip[p2$layout$name=="panel"] <- "off"

grid.arrange( p1, p2, ncol=2)


#plot(density.from, xlim = xrange, ylim = yrange, main=paste(from, '->', to), xlab = 'p', ylab = '', col='red')
#lines(density.to, col = 'blue')
#legend("topright", c(paste0('marginal: p(', from, ')'), paste0('marginal: p(', to, ')')), pch = c(1,1), col=c('red', 'blue'), border = NA)

#density.to.given.from = density(p.to.given.from)
#density.to.not.given.from = density(p.to.not.given.from)

#plot(density.to.given.from, xlim = xrange, ylim = yrange, main=paste(from, '->', to), xlab = 'p', ylab = '', col ='red')
#lines(density.to.not.given.from, col = 'blue')

dev.copy2pdf(file = 'plot.giulio.pdf')