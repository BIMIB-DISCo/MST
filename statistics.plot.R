library(lattice)
library(RColorBrewer)

performance.plot <- function(dataset,
    sample.type,
    branching) {

    if (! branching %in% c('low', 'medium', 'high', 'random_5', 'random_10', 'random_15', 'random_20')) {
        stop('branching must be "low", "medium" or "high"')
    } 
    if (! dir.exists(branching)) {
        dir.create(branching)
    }
    if (! dir.exists(paste0(branching, '/', sample.type))) {
        dir.create(paste0(branching, '/', sample.type))
    }

    sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
    sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)

    if (sample.type == "single") {
        sample = sample_sizes_single_cells
    } else if (sample.type == "multiple") {
        sample = sample_sizes_multiple_biopses
    } else {
        stop('sample.type must be "single" or "multiple"\n')
    }


    for (type in names(dataset)) {
        cat(type, ' ')
        if (! dir.exists(paste0(branching, '/', sample.type, '/', type))) {
            dir.create(paste0(branching, '/', sample.type, '/', type))
        }
        algorithms = get(type, dataset)
        for (algorithm in names(algorithms)) {
            cat('\n  ', algorithm, ' ')
            regs = get(algorithm, algorithms)
            for (reg in names(regs)) {
                cat('\n    ', reg)
                stats = get(reg, regs)

                m = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                median = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                sd = matrix(0, nrow = nrow(stats), ncol = ncol(stats))

                for (i in 1:nrow(stats)) {
                    for (j in 1:ncol(stats)) {
                        m[i,j] = stats[[i,j]]$mean
                        median[i,j] = stats[[i,j]]$median
                        sd[i,j] = stats[[i,j]]$sd
                    }
                }

                zlim = range(stats)
                zlim[1] = zlim[1] - 0.05
                zlim[2] = 1.05

                # x row noise
                # y col sample

                colnames(stats) = as.character(sample)
                rownames(stats) = as.character(1:nrow(stats))

                pdf(file=paste0(branching, '/', sample.type, '/', type, '/', algorithm, '_', reg, '.pdf'), width=8.5, height=6.5)
                print(wireframe(stats,
                    main = toupper(paste(algorithm, reg, type)),
                    drape=TRUE,
                    colorkey=TRUE,
                    xlab = 'noise',  
                    ylab = 'sample', 
                    zlab = '',
                    zlim = zlim,
                    screen = list(z = -45, x = -60),
                    scales = list(arrows = FALSE)))
                dev.off()
            }
        }
        cat('\n')
    }
}

compare.performance.plot <- function(dataset,
    sample.type,
    branching) {

    if (! branching %in% c('low', 'medium', 'high', 'random_5', 'random_10', 'random_15', 'random_20')) {
        stop('branching must be "low", "medium" or "high"')
    } 
    if (! dir.exists(branching)) {
        dir.create(branching)
    }
    if (! dir.exists(paste0(branching, '/', sample.type))) {
        dir.create(paste0(branching, '/', sample.type))
    }

    sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
    sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)

    if (sample.type == "single") {
        sample = sample_sizes_single_cells
    } else if (sample.type == "multiple") {
        sample = sample_sizes_multiple_biopses
    } else {
        stop('sample.type must be "single" or "multiple"\n')
    }

    


    for (type in names(dataset)) {
        cat(type, ' ')
        if (! dir.exists(paste0(branching, '/', sample.type, '/', type))) {
            dir.create(paste0(branching, '/', sample.type, '/', type))
        }


        results = data.frame(x = NULL, stringsAsFactors = FALSE)
        ordered.regs = NULL

        results.compare = data.frame(x = NULL, stringsAsFactors = FALSE)
        ordered.regs.compare = NULL

        results.compare.mle = data.frame(x = NULL, stringsAsFactors = FALSE)
        ordered.regs.compare.mle = NULL        

        results.compare.mltree = data.frame(x = NULL, stringsAsFactors = FALSE)
        ordered.regs.compare.mltree = NULL   

        results.compare.gabow = data.frame(x = NULL, stringsAsFactors = FALSE)
        ordered.regs.compare.gabow = NULL   

        algorithms = get(type, dataset)
        for (algorithm in names(algorithms)) {
            cat('\n  ', algorithm, ' ')
            regs = get(algorithm, algorithms)
            for (reg in names(regs)) {

                stats = get(reg, regs)

                m = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                median = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                sd = matrix(0, nrow = nrow(stats), ncol = ncol(stats))

                for (i in 1:nrow(stats)) {
                    for (j in 1:ncol(stats)) {
                        m[i,j] = stats[[i,j]]$mean
                        median[i,j] = stats[[i,j]]$median
                        sd[i,j] = stats[[i,j]]$sd
                    }
                }

                if ((algorithm == 'caprese' && reg == 'no.reg')
                    || (algorithm == 'prim' && reg == 'no.reg')
                    || (algorithm == 'edmonds')
                    || (algorithm == 'gabow')
                    || (algorithm == 'capri' && reg == 'aic')
                    || (algorithm == 'capri' && reg == 'bic')
                    || (algorithm == 'mle' && reg == 'no.reg')
                    || (algorithm == 'mltree' && reg == 'no.reg')
                    || (algorithm == 'chowliu' && reg == 'loglik')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            results = rbind(results, c(as.numeric(i),
                                as.numeric(j),
                                as.numeric(m[i,j]),
                                as.numeric(sd[i,j]),
                                as.numeric(median[i,j]),
                                paste(algorithm, reg)), stringsAsFactors = FALSE)
                        }
                    }
                    ordered.regs = c(ordered.regs, paste(algorithm, reg))
                }

                if ((algorithm == 'edmonds')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            results.compare = rbind(results.compare, c(as.numeric(i),
                                as.numeric(j),
                                as.numeric(m[i,j]),
                                as.numeric(sd[i,j]),
                                as.numeric(median[i,j]),
                                paste(algorithm, reg)), stringsAsFactors = FALSE)
                        }
                    }
                    ordered.regs.compare = c(ordered.regs.compare, paste(algorithm, reg))
                }

                if ((algorithm == 'mle' && reg == 'no.reg')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            results.compare.mle = rbind(results.compare.mle, c(as.numeric(i),
                                as.numeric(j),
                                as.numeric(m[i,j]),
                                as.numeric(sd[i,j]),
                                as.numeric(median[i,j]),
                                paste(algorithm, reg)), stringsAsFactors = FALSE)
                        }
                    }
                    ordered.regs.compare.mle = c(ordered.regs.compare.mle, paste(algorithm, reg))
                }

                if ((algorithm == 'mltree' && reg == 'no.reg')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            results.compare.mltree = rbind(results.compare.mltree, c(as.numeric(i),
                                as.numeric(j),
                                as.numeric(m[i,j]),
                                as.numeric(sd[i,j]),
                                as.numeric(median[i,j]),
                                paste(algorithm, reg)), stringsAsFactors = FALSE)
                        }
                    }
                    ordered.regs.compare.mltree = c(ordered.regs.compare.mltree, paste(algorithm, reg))
                }

                if ((algorithm == 'gabow')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            results.compare.gabow = rbind(results.compare.gabow, c(as.numeric(i),
                                as.numeric(j),
                                as.numeric(m[i,j]),
                                as.numeric(sd[i,j]),
                                as.numeric(median[i,j]),
                                paste(algorithm, reg)), stringsAsFactors = FALSE)
                        }
                    }
                    ordered.regs.compare.gabow = c(ordered.regs.compare.gabow, paste(algorithm, reg))
                }
            }
        }



        colnames(results) = c('x', 'y', 'mean', 'sd', 'median', 'algorithm')
        #colnames(results.compare) = c('x', 'y', 'mean', 'sd', 'median', 'algorithm')
        #colnames(results.compare.mle) = c('x', 'y', 'mean', 'sd', 'median', 'algorithm')
        #colnames(results.compare.mltree) = c('x', 'y', 'mean', 'sd', 'median', 'algorithm')
        colnames(results.compare.gabow) = c('x', 'y', 'mean', 'sd', 'median', 'algorithm')

        results$x = as.numeric(results$x)
        results$y = as.numeric(results$y)
        results$mean = as.numeric(results$mean)
        results$sd = as.numeric(results$sd)
        results$median = as.numeric(results$median)
        results$algorithm = as.factor(results$algorithm)

        #results.compare$x = as.numeric(results.compare$x)
        #results.compare$y = as.numeric(results.compare$y)
        #results.compare$mean = as.numeric(results.compare$mean)
        #results.compare$sd = as.numeric(results.compare$sd)
        #results.compare$median = as.numeric(results.compare$median)
        #results.compare$algorithm = as.factor(results.compare$algorithm)


        #results.compare.mle$x = as.numeric(results.compare.mle$x)
        #results.compare.mle$y = as.numeric(results.compare.mle$y)
        #results.compare.mle$mean = as.numeric(results.compare.mle$mean)
        #results.compare.mle$sd = as.numeric(results.compare.mle$sd)
        #results.compare.mle$median = as.numeric(results.compare.mle$median)
        #results.compare.mle$algorithm = as.factor(results.compare.mle$algorithm)

        #results.compare.mltree$x = as.numeric(results.compare.mltree$x)
        #results.compare.mltree$y = as.numeric(results.compare.mltree$y)
        #results.compare.mltree$mean = as.numeric(results.compare.mltree$mean)
        #results.compare.mltree$sd = as.numeric(results.compare.mltree$sd)
        #results.compare.mltree$median = as.numeric(results.compare.mltree$median)
        #results.compare.mltree$algorithm = as.factor(results.compare.mltree$algorithm)

        results.compare.gabow$x = as.numeric(results.compare.gabow$x)
        results.compare.gabow$y = as.numeric(results.compare.gabow$y)
        results.compare.gabow$mean = as.numeric(results.compare.gabow$mean)
        results.compare.gabow$sd = as.numeric(results.compare.gabow$sd)
        results.compare.gabow$median = as.numeric(results.compare.gabow$median)
        results.compare.gabow$algorithm = as.factor(results.compare.gabow$algorithm)

        print('fatto')

        # plot mix
        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'compare', '.pdf'), width=8.5, height=6.5)

        zlim = c(min(results[,'mean']) - 0.05, max(results[,'mean']) + 0.05)
        mycolors = brewer.pal(length(ordered.regs), 'Set1')
        mycolors.trans = add.alpha(mycolors, 0.7)


        print(wireframe(mean~x*y,
            data = results,
            groups = algorithm,
            col.groups = mycolors.trans,
            scales = list(arrows = FALSE,
                y = list(at = 1:length(sample), labels = as.character(sample))),
            main = toupper(paste('compare', branching, type)),
            xlab = 'noise',
            ylab = 'sample',
            zlab = '',
            screen = list(z = -45, x = -60),
            key = list(text=list(sort(ordered.regs), col = mycolors),
                    lines = list(lty = rep(1, length(ordered.regs)), col = mycolors)),
            zlim = zlim))

        dev.off()

        # plot compare

#        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'edmonds_scite_compare', '.pdf'), width=8.5, height=6.5)
#
#        zlim = c(min(results.compare[,'mean']) - 0.05, max(results.compare[,'mean']) + 0.05)
#        mycolors = brewer.pal(length(ordered.regs.compare), 'Set1')[1:length(ordered.regs.compare)]
#        mycolors.trans = add.alpha(mycolors, 0.7)
#
#        print(wireframe(mean~x*y,
#            data = results.compare,
#            groups = algorithm,
#            col.groups = mycolors.trans,
#            scales = list(arrows = FALSE,
#                y = list(at = 1:length(sample), labels = as.character(sample))),
#            main = toupper(paste('compare', branching, type)),
#            xlab = 'noise',
#            ylab = 'sample',
#            zlab = '',
#            screen = list(z = -45, x = -60),
#            key = list(text=list(sort(ordered.regs.compare), col = mycolors),
#                    lines = list(lty = rep(1, length(ordered.regs.compare)), col = mycolors)),
#            zlim = zlim))
#
#        dev.off()


        # plot compare mle

#        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'mle_scite_compare', '.pdf'), width=8.5, height=6.5)
#
#        zlim = c(min(results.compare.mle[,'mean']) - 0.05, max(results.compare.mle[,'mean']) + 0.05)
#        mycolors = brewer.pal(length(ordered.regs.compare.mle), 'Accent')[1:length(ordered.regs.compare.mle)]
#        mycolors.trans = add.alpha(mycolors, 0.7)
#
#        print(wireframe(mean~x*y,
#            data = results.compare.mle,
#            groups = algorithm,
#            col.groups = mycolors.trans,
#            scales = list(arrows = FALSE,
#                y = list(at = 1:length(sample), labels = as.character(sample))),
#            main = toupper(paste('compare', branching, type)),
#            xlab = 'noise',
#            ylab = 'sample',
#            zlab = '',
#            screen = list(z = -45, x = -60),
#            key = list(text=list(sort(ordered.regs.compare.mle), col = mycolors),
#                    lines = list(lty = rep(1, length(ordered.regs.compare.mle)), col = mycolors)),
#            zlim = zlim))
#
#        dev.off()
#
#        # plot compare mltree
#
#        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'mltree_scite_compare', '.pdf'), width=8.5, height=6.5)
#
#        zlim = c(min(results.compare.mltree[,'mean']) - 0.05, max(results.compare.mltree[,'mean']) + 0.05)
#        mycolors = brewer.pal(length(ordered.regs.compare.mltree), 'Accent')[1:length(ordered.regs.compare.mltree)]
#        mycolors.trans = add.alpha(mycolors, 0.7)
#
#        print(wireframe(mean~x*y,
#            data = results.compare.mltree,
#            groups = algorithm,
#            col.groups = mycolors.trans,
#            scales = list(arrows = FALSE,
#                y = list(at = 1:length(sample), labels = as.character(sample))),
#            main = toupper(paste('compare', branching, type)),
#            xlab = 'noise',
#            ylab = 'sample',
#            zlab = '',
#            screen = list(z = -45, x = -60),
#            key = list(text=list(sort(ordered.regs.compare.mltree), col = mycolors),
#                    lines = list(lty = rep(1, length(ordered.regs.compare.mltree)), col = mycolors)),
#            zlim = zlim))
#
#        dev.off()

        # plot compare gabow

        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'gabow_scite_compare', '.pdf'), width=8.5, height=6.5)

        zlim = c(min(results.compare.gabow[,'mean']) - 0.05, max(results.compare.gabow[,'mean']) + 0.05)
        mycolors = brewer.pal(length(ordered.regs.compare.gabow), 'Set1')[1:length(ordered.regs.compare.gabow)]
        mycolors.trans = add.alpha(mycolors, 0.7)

        print(results.compare.gabow)

        print(wireframe(mean~x*y,
            data = results.compare.gabow,
            groups = algorithm,
            col.groups = mycolors.trans,
            scales = list(arrows = FALSE,
                y = list(at = 1:length(sample), labels = as.character(sample))),
            main = toupper(paste('compare', branching, type)),
            xlab = 'noise',
            ylab = 'sample',
            zlab = '',
            screen = list(z = -45, x = -60),
            key = list(text=list(sort(ordered.regs.compare.gabow), col = mycolors),
                    lines = list(lty = rep(1, length(ordered.regs.compare.gabow)), col = mycolors)),
            zlim = zlim))

        dev.off()

        # plot compare median

#        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'edmonds_scite_compare_median', '.pdf'), width=8.5, height=6.5)
#
#        zlim = c(min(results.compare[,'median']) - 0.05, max(results.compare[,'median']) + 0.05)
#        mycolors = brewer.pal(length(ordered.regs.compare), 'Set1')[1:length(ordered.regs.compare)]
#        mycolors.trans = add.alpha(mycolors, 0.7)
#
#        print(wireframe(median~x*y,
#            data = results.compare,
#            groups = algorithm,
#            col.groups = mycolors.trans,
#            scales = list(arrows = FALSE,
#                y = list(at = 1:length(sample), labels = as.character(sample))),
#            main = toupper(paste('compare', branching, type)),
#            xlab = 'noise',
#            ylab = 'sample',
#            zlab = '',
#            screen = list(z = -45, x = -60),
#            key = list(text=list(sort(ordered.regs.compare), col = mycolors),
#                    lines = list(lty = rep(1, length(ordered.regs.compare)), col = mycolors)),
#            zlim = zlim))
#
#        dev.off()


#        # plot compare mle median
#
#        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'mle_scite_compare_median', '.pdf'), width=8.5, height=6.5)
#
#        zlim = c(min(results.compare.mle[,'median']) - 0.05, max(results.compare.mle[,'median']) + 0.05)
#        mycolors = brewer.pal(length(ordered.regs.compare.mle), 'Accent')[1:length(ordered.regs.compare.mle)]
#        mycolors.trans = add.alpha(mycolors, 0.7)
#
#        print(wireframe(median~x*y,
#            data = results.compare.mle,
#            groups = algorithm,
#            col.groups = mycolors.trans,
#            scales = list(arrows = FALSE,
#                y = list(at = 1:length(sample), labels = as.character(sample))),
#            main = toupper(paste('compare', branching, type)),
#            xlab = 'noise',
#            ylab = 'sample',
#            zlab = '',
#            screen = list(z = -45, x = -60),
#            key = list(text=list(sort(ordered.regs.compare.mle), col = mycolors),
#                    lines = list(lty = rep(1, length(ordered.regs.compare.mle)), col = mycolors)),
#            zlim = zlim))
#
#        dev.off()
#
#        # plot compare mltree median
#
#        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'mltree_scite_compare_median', '.pdf'), width=8.5, height=6.5)
#
#        zlim = c(min(results.compare.mltree[,'median']) - 0.05, max(results.compare.mltree[,'median']) + 0.05)
#        mycolors = brewer.pal(length(ordered.regs.compare.mltree), 'Accent')[1:length(ordered.regs.compare.mltree)]
#        mycolors.trans = add.alpha(mycolors, 0.7)
#
#        print(wireframe(median~x*y,
#            data = results.compare.mltree,
#            groups = algorithm,
#            col.groups = mycolors.trans,
#            scales = list(arrows = FALSE,
#                y = list(at = 1:length(sample), labels = as.character(sample))),
#            main = toupper(paste('compare', branching, type)),
#            xlab = 'noise',
#            ylab = 'sample',
#            zlab = '',
#            screen = list(z = -45, x = -60),
#            key = list(text=list(sort(ordered.regs.compare.mltree), col = mycolors),
#                    lines = list(lty = rep(1, length(ordered.regs.compare.mltree)), col = mycolors)),
#            zlim = zlim))
#
#        dev.off()

        # plot compare gabow median

        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'gabow_scite_compare_median', '.pdf'), width=8.5, height=6.5)

        zlim = c(min(results.compare.gabow[,'median']) - 0.05, max(results.compare.gabow[,'median']) + 0.05)
        mycolors = brewer.pal(length(ordered.regs.compare.gabow), 'Set1')[1:length(ordered.regs.compare.gabow)]
        mycolors.trans = add.alpha(mycolors, 0.7)

        print(wireframe(median~x*y,
            data = results.compare.gabow,
            groups = algorithm,
            col.groups = mycolors.trans,
            scales = list(arrows = FALSE,
                y = list(at = 1:length(sample), labels = as.character(sample))),
            main = toupper(paste('compare', branching, type)),
            xlab = 'noise',
            ylab = 'sample',
            zlab = '',
            screen = list(z = -45, x = -60),
            key = list(text=list(sort(ordered.regs.compare.gabow), col = mycolors),
                    lines = list(lty = rep(1, length(ordered.regs.compare.gabow)), col = mycolors)),
            zlim = zlim))

        dev.off()

        cat('\n')
    }
}


compare.performance.plot.2d <- function(dataset,
    sample.type,
    branching) {

    if (! branching %in% c('low', 'medium', 'high', 'random_5', 'random_10', 'random_15', 'random_20')) {
        stop('branching must be "low", "medium" or "high"')
    } 
    if (! dir.exists(branching)) {
        dir.create(branching)
    }
    if (! dir.exists(paste0(branching, '/', sample.type))) {
        dir.create(paste0(branching, '/', sample.type))
    }

    sample_sizes_single_cells = c(10, 50, 100, 250)
    sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)

    if (sample.type == "single") {
        sample = sample_sizes_single_cells
        select.sample = c(10, 50, 100, 250)
    } else if (sample.type == "multiple") {
        sample = sample_sizes_multiple_biopses
        select.sample = c(5, 20, 100)
    } else {
        stop('sample.type must be "single" or "multiple"\n')
    }

    


    for (type in names(dataset)) {
        cat(type, ' ')
        if (! dir.exists(paste0(branching, '/', sample.type, '/', type))) {
            dir.create(paste0(branching, '/', sample.type, '/', type))
        }


        results.c1 = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.c2 = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.c3 = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.c4 = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.c5 = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.c6 = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.c7 = data.frame(x = NULL, stringsAsFactors = FALSE)
        results.c8 = data.frame(x = NULL, stringsAsFactors = FALSE)
        ordered.regs.c1 = NULL
        ordered.regs.c2 = NULL
        ordered.regs.c3 = NULL
        ordered.regs.c4 = NULL
        ordered.regs.c5 = NULL
        ordered.regs.c6 = NULL
        ordered.regs.c7 = NULL
        ordered.regs.c8 = NULL

        algorithms = get(type, dataset)
        for (algorithm in names(algorithms)) {
            cat('\n  ', algorithm, ' ')
            regs = get(algorithm, algorithms)
            for (reg in names(regs)) {

                stats = get(reg, regs)

                m = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                median = matrix(0, nrow = nrow(stats), ncol = ncol(stats))
                sd = matrix(0, nrow = nrow(stats), ncol = ncol(stats))

                for (i in 1:nrow(stats)) {
                    for (j in 1:ncol(stats)) {
                        m[i,j] = stats[[i,j]]$mean
                        median[i,j] = stats[[i,j]]$median
                        sd[i,j] = stats[[i,j]]$sd
                    }
                }

                # cfg 1 - mix
                if ((algorithm == 'caprese' && reg == 'no.reg')
                    || (algorithm == 'prim' && reg == 'no.reg')
                    || (algorithm == 'edmonds')
                    || (algorithm == 'gabow')
                    || (algorithm == 'gabow_no_rising')
                    || (algorithm == 'chowliu' && reg == 'loglik')
                    || (algorithm == 'capri' && reg == 'bic')
                    || (algorithm == 'capri' && reg == 'aic')
                    || (algorithm == 'mle' && reg == 'no.reg')
                    || (algorithm == 'mltree' && reg == 'no.reg')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c1 = rbind(results.c1, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    paste(algorithm, reg)), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c1 = c(ordered.regs.c1, paste(algorithm, reg))
                }

                # cfg 2 - all prim
                if (algorithm == 'prim') {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c2 = rbind(results.c2, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    reg), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c2 = c(ordered.regs.c2, reg)
                }

                # cfg 3 - all edmonds
                if (algorithm == 'edmonds' && !reg %in% c('bic', 'aic', 'or')) {
                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c3 = rbind(results.c3, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    reg), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c3 = c(ordered.regs.c3, reg)
                }

                # cfg 4 - all chowliu
                if (algorithm == 'chowliu') {
                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c4 = rbind(results.c4, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    reg), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c4 = c(ordered.regs.c4, reg)
                }

                # cfg 5 - edmonds vs scite
                if ((algorithm == 'edmonds')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c5 = rbind(results.c5, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    paste(algorithm, reg)), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c5 = c(ordered.regs.c5, paste(algorithm, reg))
                }


                # cfg 6 - mle vs scite
                if ((algorithm == 'mle' && !reg %in% c('bic', 'aic', 'cmi2', 'cmi3', 'or'))
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c6 = rbind(results.c6, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    paste(algorithm, reg)), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c6 = c(ordered.regs.c6, paste(algorithm, reg))
                }

                # cfg 2 - mltree vs scite
                if ((algorithm == 'mltree' && !reg %in% c('bic', 'aic', 'cmi2', 'cmi3', 'or'))
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c2 = rbind(results.c2, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    paste(algorithm, reg)), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c2 = c(ordered.regs.c2, paste(algorithm, reg))
                }

                # cfg 7 - gabow vs scite
                if ((algorithm == 'gabow')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c7 = rbind(results.c7, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    paste(algorithm, reg)), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c7 = c(ordered.regs.c7, paste(algorithm, reg))
                }

                # cfg 8 - gabow no rais vs scite
                if ((algorithm == 'gabow_no_rising')
                    || (algorithm == 'scite' && reg == 'no.reg')) {

                    cat('\n    ', reg)
                    for (i in 1:nrow(m)) {
                        for(j in 1:ncol(m)) {
                            if (sample[j] %in% select.sample) {
                                results.c8 = rbind(results.c8, c(i,
                                    m[i,j],
                                    sd[i,j],
                                    median[i,j],
                                    sample[j],
                                    paste(algorithm, reg)), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c8 = c(ordered.regs.c8, paste(algorithm, reg))
                }

            }
        }

        colnames(results.c1) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        #colnames(results.c2) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        #colnames(results.c3) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        #colnames(results.c5) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        #colnames(results.c6) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        colnames(results.c7) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        colnames(results.c8) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')

        results.c1$x = as.numeric(results.c1$x)
        results.c1$mean = as.numeric(results.c1$mean)
        results.c1$sd = as.numeric(results.c1$sd)
        results.c1$median = as.numeric(results.c1$median)
        results.c1$sample = as.numeric(results.c1$sample)
        results.c1$algorithm = as.factor(results.c1$algorithm)

        #results.c2$x = as.numeric(results.c2$x)
        #results.c2$mean = as.numeric(results.c2$mean)
        #results.c2$sd = as.numeric(results.c2$sd)
        #results.c2$median = as.numeric(results.c2$median)
        #results.c2$sample = as.numeric(results.c2$sample)
        #results.c2$algorithm = as.factor(results.c2$algorithm)

        #results.c3$x = as.numeric(results.c3$x)
        #results.c3$mean = as.numeric(results.c3$mean)
        #results.c3$sd = as.numeric(results.c3$sd)
        #results.c3$median = as.numeric(results.c3$median)
        #results.c3$sample = as.numeric(results.c3$sample)
        #results.c3$algorithm = as.factor(results.c3$algorithm)

        #results.c5$x = as.numeric(results.c5$x)
        #results.c5$mean = as.numeric(results.c5$mean)
        #results.c5$sd = as.numeric(results.c5$sd)
        #results.c5$median = as.numeric(results.c5$median)
        #results.c5$sample = as.numeric(results.c5$sample)
        #results.c5$algorithm = as.factor(results.c5$algorithm)

        #results.c6$x = as.numeric(results.c6$x)
        #results.c6$mean = as.numeric(results.c6$mean)
        #results.c6$sd = as.numeric(results.c6$sd)
        #results.c6$median = as.numeric(results.c6$median)
        #results.c6$sample = as.numeric(results.c6$sample)
        #results.c6$algorithm = as.factor(results.c6$algorithm)

        results.c7$x = as.numeric(results.c7$x)
        results.c7$mean = as.numeric(results.c7$mean)
        results.c7$sd = as.numeric(results.c7$sd)
        results.c7$median = as.numeric(results.c7$median)
        results.c7$sample = as.numeric(results.c7$sample)
        results.c7$algorithm = as.factor(results.c7$algorithm)

        results.c8$x = as.numeric(results.c8$x)
        results.c8$mean = as.numeric(results.c8$mean)
        results.c8$sd = as.numeric(results.c8$sd)
        results.c8$median = as.numeric(results.c8$median)
        results.c8$sample = as.numeric(results.c8$sample)
        results.c8$algorithm = as.factor(results.c8$algorithm)


        cat('\ndone\n')

        for (s in select.sample) {

            cat(s, '\n')

            # mix res 1
            cat('\tmix res mean\n')

            res = results.c1[which(results.c1$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mix', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'mean']) - 0.05, max(res[,'mean']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c1), 'Paired')
            print(xyplot(mean~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('compare sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c1), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c1)), col = mycolors)),
                ylim = ylim))

            dev.off()

            # mix res median 1
            cat('\tmix res median\n')

            res = results.c1[which(results.c1$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mix_median', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c1), 'Paired')
            print(xyplot(median~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('compare sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c1), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c1)), col = mycolors)),
                ylim = ylim))

            dev.off()

            # prim res 2

            #res = results.c2[which(results.c2$sample == s), ]
            #pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'prim', '.pdf'), width=8.5, height=6.5)
            #ylim = c(min(res[,'mean']) - 0.05, max(res[,'mean']) + 0.05)
            #mycolors = brewer.pal(length(ordered.regs.c2), 'Set1')[1:length(ordered.regs.c2)]
            #print(xyplot(mean~x,
            #    grid = TRUE,
            #    data = res,
            #    groups = algorithm,
            #    col = mycolors,
            #    main = toupper(paste('prim sample:', s, branching, type)),
            #    xlab = 'noise',
            #    ylab = '',
            #    type = c("o"),
            #    key = list(text=list(sort(ordered.regs.c2), col = mycolors),
            #            lines = list(lty = rep(1, length(ordered.regs.c2)), col = mycolors)),
            #    ylim = ylim))
            #dev.off()


            # edmonds res 3

            #res = results.c3[which(results.c3$sample == s), ]
            #pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'edmonds', '.pdf'), width=8.5, height=6.5)
            #ylim = c(min(res[,'mean']) - 0.05, max(res[,'mean']) + 0.05)
            #mycolors = brewer.pal(length(ordered.regs.c3), 'Set1')[1:length(ordered.regs.c3)]
            #print(xyplot(mean~x,
            #    grid = TRUE,
            #    data = res,
            #    groups = algorithm,
            #    col = mycolors,
            #    main = toupper(paste('edmonds sample:', s, branching, type)),
            #    xlab = 'noise',
            #    ylab = '',
            #    type = c("o"),
            #    key = list(text=list(sort(ordered.regs.c3), col = mycolors),
            #            lines = list(lty = rep(1, length(ordered.regs.c3)), col = mycolors)),
            #    ylim = ylim))
            #dev.off()

            # chowliu res 4

            #res = results.c4[which(results.c4$sample == s), ]
            #pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'chowliu', '.pdf'), width=8.5, height=6.5)
            #ylim = c(min(res[,'mean']) - 0.05, max(res[,'mean']) + 0.05)
            #mycolors = brewer.pal(length(ordered.regs.c4), 'Set1')[1:length(ordered.regs.c4)]
            #print(xyplot(mean~x,
            #    grid = TRUE,
            #    data = res,
            #    groups = algorithm,
            #    col = mycolors,
            #    main = toupper(paste('chowliu sample:', s, branching, type)),
            #    xlab = 'noise',
            #    ylab = '',
            #    type = c("o"),
            #    key = list(text=list(sort(ordered.regs.c4), col = mycolors),
            #            lines = list(lty = rep(1, length(ordered.regs.c4)), col = mycolors)),
            #    ylim = ylim))
            #dev.off()

            # edmonds vs scite 5


            prepanel.ci <- function(x, y, ly, uy, subscripts, ...) {
                x = as.numeric(x)
                ly = as.numeric(ly[subscripts])
                uy = as.numeric(uy[subscripts])
                list(ylim = range(y, uy, ly, finite = TRUE))
            }

            panel.ci <- function(x, y, ly, uy, subscripts, pch = 16, col = col, ...){
                x = as.numeric(x)
                y = as.numeric(y)
                ly = as.numeric(ly[subscripts])
                uy = as.numeric(uy[subscripts])
                panel.arrows(x, ly, x, uy, col = col,
                             length = 0.1, unit = "native",
                             angle = 90, code = 3)
                panel.xyplot(x, y, pch = pch, ...)
            }

            # edmonds vs scite

#            res = results.c5[which(results.c5$sample == s), ]
#            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'edmods_scite', '.pdf'), width=8.5, height=6.5)
#            ylim = c(min(res[,'mean']) - max(res[,'sd']) - 0.05, max(res[,'mean']) + max(res[,'sd']) + 0.05)
#            mycolors = brewer.pal(length(ordered.regs.c5), 'Set1')[1:length(ordered.regs.c5)]
#            print(xyplot(mean~x,
#                grid = TRUE,
#                data = res,
#                groups = algorithm,
#                ly = res$mean - res$sd,
#                uy = res$mean + res$sd,
#                prepanel = prepanel.ci,
#                panel = panel.superpose,
#                panel.groups = panel.ci,
#                col = mycolors,
#                main = toupper(paste('mean edmonds vs scite sample:', s, branching, type)),
#                xlab = 'noise',
#                ylab = '',
#                type = c("o"),
#                key = list(text=list(sort(ordered.regs.c5), col = mycolors),
#                        lines = list(lty = rep(1, length(ordered.regs.c5)), col = mycolors)),
#                ylim = ylim))
#            #ggplot(res, aes(x=x, y=mean, colour=algorithm))
#            #stop()
#            dev.off()


            # mle vs scite

 #           res = results.c6[which(results.c6$sample == s), ]
 #           pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mle_scite', '.pdf'), width=8.5, height=6.5)
 #           ylim = c(min(res[,'mean']) - max(res[,'sd']) - 0.05, max(res[,'mean']) + max(res[,'sd']) + 0.05)
 #           mycolors = brewer.pal(length(ordered.regs.c6), 'Set1')[1:length(ordered.regs.c6)]
 #           print(xyplot(mean~x,
 #               grid = TRUE,
 #               data = res,
 #               groups = algorithm,
 #               ly = res$mean - res$sd,
 #               uy = res$mean + res$sd,
 #               prepanel = prepanel.ci,
 #               panel = panel.superpose,
 #               panel.groups = panel.ci,
 #               col = mycolors,
 #               main = toupper(paste('mean mle vs scite sample:', s, branching, type)),
 #               xlab = 'noise',
 #               ylab = '',
 #               type = c("o"),
 #               key = list(text=list(sort(ordered.regs.c6), col = mycolors),
 #                       lines = list(lty = rep(1, length(ordered.regs.c6)), col = mycolors)),
 #               ylim = ylim))
 #           #ggplot(res, aes(x=x, y=mean, colour=algorithm))
 #           #stop()
 #           dev.off()

            # mltree vs scite

#            res = results.c2[which(results.c2$sample == s), ]
#            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mltree_scite', '.pdf'), width=8.5, height=6.5)
#            ylim = c(min(res[,'mean']) - max(res[,'sd']) - 0.05, max(res[,'mean']) + max(res[,'sd']) + 0.05)
#            mycolors = brewer.pal(length(ordered.regs.c2), 'Set1')[1:length(ordered.regs.c2)]
#            print(xyplot(mean~x,
#                grid = TRUE,
#                data = res,
#                groups = algorithm,
#                ly = res$mean - res$sd,
#                uy = res$mean + res$sd,
#                prepanel = prepanel.ci,
#                panel = panel.superpose,
#                panel.groups = panel.ci,
#                col = mycolors,
#                main = toupper(paste('mean mltree vs scite sample:', s, branching, type)),
#                xlab = 'noise',
#                ylab = '',
#                type = c("o"),
#                key = list(text=list(sort(ordered.regs.c2), col = mycolors),
#                        lines = list(lty = rep(1, length(ordered.regs.c2)), col = mycolors)),
#                ylim = ylim))
#            #ggplot(res, aes(x=x, y=mean, colour=algorithm))
#            #stop()
#            dev.off()


            # gabow vs scite
            cat('\tgabow vs scite\n')

            res = results.c7[which(results.c7$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'gabow_scite', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'mean']) - max(res[,'sd']) - 0.05, max(res[,'mean']) + max(res[,'sd']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c7), 'Paired')[1:length(ordered.regs.c7)]
            print(xyplot(mean~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                ly = res$mean - res$sd,
                uy = res$mean + res$sd,
                prepanel = prepanel.ci,
                panel = panel.superpose,
                panel.groups = panel.ci,
                col = mycolors,
                main = toupper(paste('mean gabow vs scite sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c7), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c7)), col = mycolors)),
                ylim = ylim))
            #ggplot(res, aes(x=x, y=mean, colour=algorithm))
            #stop()
            dev.off()


            # gabow f vs scite
            cat('\tgabow f vs scite\n')

            res = results.c8[which(results.c8$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'gabow_no_raising_scite', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'mean']) - max(res[,'sd']) - 0.05, max(res[,'mean']) + max(res[,'sd']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c8), 'Paired')[1:length(ordered.regs.c8)]
            print(xyplot(mean~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                ly = res$mean - res$sd,
                uy = res$mean + res$sd,
                prepanel = prepanel.ci,
                panel = panel.superpose,
                panel.groups = panel.ci,
                col = mycolors,
                main = toupper(paste('mean gabow no raising vs scite sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c8), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c8)), col = mycolors)),
                ylim = ylim))
            #ggplot(res, aes(x=x, y=mean, colour=algorithm))
            #stop()
            dev.off()

            # mltree vs scite

#            res = results.c2[which(results.c2$sample == s), ]
#            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mltree_scite', '.pdf'), width=8.5, height=6.5)
#            ylim = c(min(res[,'mean']) - max(res[,'sd']) - 0.05, max(res[,'mean']) + max(res[,'sd']) + 0.05)
#            mycolors = brewer.pal(length(ordered.regs.c2), 'Set1')[1:length(ordered.regs.c2)]
#            print(xyplot(mean~x,
#                grid = TRUE,
#                data = res,
#                groups = algorithm,
#                ly = res$mean - res$sd,
#                uy = res$mean + res$sd,
#                prepanel = prepanel.ci,
#                panel = panel.superpose,
#                panel.groups = panel.ci,
#                col = mycolors,
#                main = toupper(paste('mean mltree vs scite sample:', s, branching, type)),
#                xlab = 'noise',
#                ylab = '',
#                type = c("o"),
#                key = list(text=list(sort(ordered.regs.c2), col = mycolors),
#                        lines = list(lty = rep(1, length(ordered.regs.c2)), col = mycolors)),
#                ylim = ylim))
#            #ggplot(res, aes(x=x, y=mean, colour=algorithm))
#            #stop()
#            dev.off()
            # edmonds vs scite 5 median

#            res = results.c5[which(results.c5$sample == s), ]
#            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'edmods_scite_median', '.pdf'), width=8.5, height=6.5)
#            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
#            mycolors = brewer.pal(length(ordered.regs.c5), 'Set1')[1:length(ordered.regs.c5)]
#            print(xyplot(median~x,
#                grid = TRUE,
#                data = res,
#                groups = algorithm,
#                col = mycolors,
#                main = toupper(paste('median edmonds vs scite sample:', s, branching, type)),
#                xlab = 'noise',
#                ylab = '',
#                type = c("o"),
#                key = list(text=list(sort(ordered.regs.c5), col = mycolors),
#                        lines = list(lty = rep(1, length(ordered.regs.c5)), col = mycolors)),
#                ylim = ylim))
#            dev.off()


#            # mle vs scite 5 median
#
#            res = results.c6[which(results.c6$sample == s), ]
#            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mle_scite_median', '.pdf'), width=8.5, height=6.5)
#            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
#            mycolors = brewer.pal(length(ordered.regs.c6), 'Set1')[1:length(ordered.regs.c6)]
#            print(xyplot(median~x,
#                grid = TRUE,
#                data = res,
#                groups = algorithm,
#                col = mycolors,
#                main = toupper(paste('median mle vs scite sample:', s, branching, type)),
#                xlab = 'noise',
#                ylab = '',
#                type = c("o"),
#                key = list(text=list(sort(ordered.regs.c6), col = mycolors),
#                        lines = list(lty = rep(1, length(ordered.regs.c6)), col = mycolors)),
#                ylim = ylim))
#            dev.off()
#
#             # mltree vs scite 5 median
#
#            res = results.c2[which(results.c2$sample == s), ]
#            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mltree_scite_median', '.pdf'), width=8.5, height=6.5)
#            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
#            mycolors = brewer.pal(length(ordered.regs.c2), 'Set1')[1:length(ordered.regs.c2)]
#            print(xyplot(median~x,
#                grid = TRUE,
#                data = res,
#                groups = algorithm,
#                col = mycolors,
#                main = toupper(paste('median mltree vs scite sample:', s, branching, type)),
#                xlab = 'noise',
#                ylab = '',
#                type = c("o"),
#                key = list(text=list(sort(ordered.regs.c2), col = mycolors),
#                        lines = list(lty = rep(1, length(ordered.regs.c2)), col = mycolors)),
#                ylim = ylim))
#            dev.off()


            # gabow vs scite 5 median
            cat('\tgabow vs scite median\n')

            res = results.c7[which(results.c7$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'gabow_scite_median', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c7), 'Paired')[1:length(ordered.regs.c7)]
            print(xyplot(median~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('median gabow vs scite sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c7), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c7)), col = mycolors)),
                ylim = ylim))
            dev.off()

            # gabow f vs scite 5 median
            cat('\tgabow f vs scite median\n')

            res = results.c8[which(results.c8$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'gabow_no_rais_scite_median', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c8), 'Paired')[1:length(ordered.regs.c8)]
            print(xyplot(median~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('median gabow vs scite sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c8), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c8)), col = mycolors)),
                ylim = ylim))
            dev.off()


        }
        cat('\n')
    }
}




add.alpha <- function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    r = apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))
    return(as.vector(r))  
}
