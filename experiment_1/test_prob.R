##############################################################################
###
### MST
###
### Test Probability
###
##############################################################################
### Copyright (c) 2015-2018, The TRONCO Team (www.troncopackage.org)
### email: tronco@disco.unimib.it
### All rights reserved. This program and the accompanying materials
### are made available under the terms of the GNU GPL v3.0
### which accompanies this distribution
##############################################################################


library(TRONCO)


ratio_left = c()
ratio_right = c()


for (i in 1:100) {
    exp = dataset.single.cells.convergent[[1, i]][[1]]
    prob = exp$probabilities
    p_root = exp$root_prob
    dataset = exp$dataset
    epos = exp$epos
    eneg = exp$eneg
    adj.matrix = array(1, c(ncol(dataset), ncol(dataset)))
    colnames(adj.matrix) = colnames(dataset)
    rownames(adj.matrix) = colnames(dataset)
    scores = TRONCO:::get.dag.scores( dataset, adj.matrix, 0, 0)
    marginal = scores$marginal.probs
    joint = scores$joint.probs

    left = joint[[2,3]] / (marginal[[2]] * marginal[[3]])
    ratio_left = c(ratio_left, left)

    right = joint[[3,4]] / (marginal[[3]] * marginal[[4]])
    ratio_right = c(ratio_right, right)
}

results = data.frame(x = NULL, stringsAsFactors = FALSE)

for (val in ratio_left) {
    results = rbind(results, c(val, 'left'), stringsAsFactors = FALSE)
}

for (val in ratio_right) {
    results = rbind(results, c(val, 'right'), stringsAsFactors = FALSE)
}

colnames(results) = c('value', 'direction')
results$direction = as.factor(results$direction)
results$value = as.numeric(results$value)

experiment.palette = c('left' = 'red', 'right' = 'blue')

experiment.names = c('left' = '2->5 3->5', 'right' = '3->6 4->6')


p = ggplot(results, aes(direction, value)) + 
    #p = ggplot(res.values, aes(x = reorder(noise, value, FUN = median), y = value, fill = algorithm)) + 
    #scale_fill_manual(values = experiment.palette, labels = experiment.names) +
    scale_x_discrete(label=experiment.names) +
    #expand_limits(y=c(0.5,1)) +
    theme(axis.text.y = element_text(size = 12), legend.position = 'none') +
    #scale_y_continuous(label = function(x){paste0(x, ' (', as.integer(x/connections*100), '%)')}) +
    ylab("P(A,B) / (P(A) * P(B))") +
    xlab('???') +
    #geom_boxplot(outlier.size = 0) +
    geom_jitter(aes(colour = direction), position = position_jitterdodge(jitter.height = 0, jitter.width = .5)) +
    ggtitle('???')

print(p)
dev.copy2pdf(file = 'test_prob.pdf')

### end of file -- test_prob.R
