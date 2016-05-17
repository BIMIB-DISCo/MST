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


                print(stats)

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
                    || (algorithm == 'edmonds' && reg == 'no.reg')
                    || (algorithm == 'capri' && reg == 'aic')
                    || (algorithm == 'capri' && reg == 'bic')
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

                if ((algorithm == 'edmonds' && reg == 'no.reg')
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
            }
        }


        print('done!!!')

        colnames(results) = c('x', 'y', 'mean', 'sd', 'median', 'algorithm')
        colnames(results.compare) = c('x', 'y', 'mean', 'sd', 'median', 'algorithm')

        results$x = as.numeric(results$x)
        results$y = as.numeric(results$y)
        results$mean = as.numeric(results$mean)
        results$sd = as.numeric(results$sd)
        results$median = as.numeric(results$median)
        results$algorithm = as.factor(results$algorithm)

        results.compare$x = as.numeric(results.compare$x)
        results.compare$y = as.numeric(results.compare$y)
        results.compare$mean = as.numeric(results.compare$mean)
        results.compare$sd = as.numeric(results.compare$sd)
        results.compare$median = as.numeric(results.compare$median)
        results.compare$algorithm = as.factor(results.compare$algorithm)

        print('ok')

        # plot mix
        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'compare', '.pdf'), width=8.5, height=6.5)

        zlim = c(min(results[,'mean']) - 0.05, max(results[,'mean']) + 0.05)
        mycolors = brewer.pal(length(ordered.regs), 'Accent')
        mycolors.trans = add.alpha(mycolors, 0.7)

        print('asd')

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

        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'edmonds_scite_compare', '.pdf'), width=8.5, height=6.5)

        zlim = c(min(results.compare[,'mean']) - 0.05, max(results.compare[,'median']) + 0.05)
        mycolors = brewer.pal(length(ordered.regs.compare), 'Accent')[1:length(ordered.regs.compare)]
        mycolors.trans = add.alpha(mycolors, 0.7)

        print(wireframe(mean~x*y,
            data = results.compare,
            groups = algorithm,
            col.groups = mycolors.trans,
            scales = list(arrows = FALSE,
                y = list(at = 1:length(sample), labels = as.character(sample))),
            main = toupper(paste('compare', branching, type)),
            xlab = 'noise',
            ylab = 'sample',
            zlab = '',
            screen = list(z = -45, x = -60),
            key = list(text=list(sort(ordered.regs.compare), col = mycolors),
                    lines = list(lty = rep(1, length(ordered.regs.compare)), col = mycolors)),
            zlim = zlim))

        dev.off()

        # plot compare median

        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'edmonds_scite_compare_median', '.pdf'), width=8.5, height=6.5)

        zlim = c(min(results.compare[,'median']) - 0.05, max(results.compare[,'median']) + 0.05)
        mycolors = brewer.pal(length(ordered.regs.compare), 'Accent')[1:length(ordered.regs.compare)]
        mycolors.trans = add.alpha(mycolors, 0.7)

        print(wireframe(median~x*y,
            data = results.compare,
            groups = algorithm,
            col.groups = mycolors.trans,
            scales = list(arrows = FALSE,
                y = list(at = 1:length(sample), labels = as.character(sample))),
            main = toupper(paste('compare', branching, type)),
            xlab = 'noise',
            ylab = 'sample',
            zlab = '',
            screen = list(z = -45, x = -60),
            key = list(text=list(sort(ordered.regs.compare), col = mycolors),
                    lines = list(lty = rep(1, length(ordered.regs.compare)), col = mycolors)),
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

    sample_sizes_single_cells = c(10, 25, 50, 75, 100, 150, 200, 250, 500, 1000)
    sample_sizes_multiple_biopses = c(5, 6, 7, 8, 9, 10, 15, 20, 50, 100)

    if (sample.type == "single") {
        sample = sample_sizes_single_cells
        select.sample = c(10, 100, 1000)
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
        ordered.regs.c1 = NULL
        ordered.regs.c2 = NULL
        ordered.regs.c3 = NULL
        ordered.regs.c4 = NULL
        ordered.regs.c5 = NULL

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
                    || (algorithm == 'edmonds' && reg == 'no.reg')
                    || (algorithm == 'chowliu' && reg == 'loglik')
                    || (algorithm == 'capri' && reg == 'bic')
                    || (algorithm == 'capri' && reg == 'aic')
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
                if ((algorithm == 'edmonds' && !reg %in% c('bic', 'aic', 'cmi2', 'cmi3', 'or'))
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

            }
        }

        colnames(results.c1) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        colnames(results.c2) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        colnames(results.c3) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        colnames(results.c4) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')
        colnames(results.c5) = c('x', 'mean', 'sd', 'median', 'sample', 'algorithm')

        results.c1$x = as.numeric(results.c1$x)
        results.c1$mean = as.numeric(results.c1$mean)
        results.c1$sd = as.numeric(results.c1$sd)
        results.c1$median = as.numeric(results.c1$median)
        results.c1$sample = as.numeric(results.c1$sample)
        results.c1$algorithm = as.factor(results.c1$algorithm)

        results.c2$x = as.numeric(results.c2$x)
        results.c2$mean = as.numeric(results.c2$mean)
        results.c2$sd = as.numeric(results.c2$sd)
        results.c2$median = as.numeric(results.c2$median)
        results.c2$sample = as.numeric(results.c2$sample)
        results.c2$algorithm = as.factor(results.c2$algorithm)

        results.c3$x = as.numeric(results.c3$x)
        results.c3$mean = as.numeric(results.c3$mean)
        results.c3$sd = as.numeric(results.c3$sd)
        results.c3$median = as.numeric(results.c3$median)
        results.c3$sample = as.numeric(results.c3$sample)
        results.c3$algorithm = as.factor(results.c3$algorithm)

        results.c4$x = as.numeric(results.c4$x)
        results.c4$mean = as.numeric(results.c4$mean)
        results.c4$sd = as.numeric(results.c4$sd)
        results.c4$median = as.numeric(results.c4$median)
        results.c4$sample = as.numeric(results.c4$sample)
        results.c4$algorithm = as.factor(results.c4$algorithm)

        results.c5$x = as.numeric(results.c5$x)
        results.c5$mean = as.numeric(results.c5$mean)
        results.c5$sd = as.numeric(results.c5$sd)
        results.c5$median = as.numeric(results.c5$median)
        results.c5$sample = as.numeric(results.c5$sample)
        results.c5$algorithm = as.factor(results.c5$algorithm)


        for (s in select.sample) {

            # mix res 1

            res = results.c1[which(results.c1$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mix', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'mean']) - 0.05, max(res[,'mean']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c1), 'Dark2')
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

            res = results.c1[which(results.c1$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mix_median', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c1), 'Dark2')
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

            res = results.c5[which(results.c5$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'edmods_scite', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'mean']) - max(res[,'sd']) - 0.05, max(res[,'mean']) + max(res[,'sd']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c5), 'Set1')[1:length(ordered.regs.c5)]
            #print(xyplot(mean~x,
            #    grid = TRUE,
            #    data = res,
            #    groups = algorithm,
            #    col = mycolors,
            #    main = toupper(paste('mean edmonds vs scite sample:', s, branching, type)),
            #    xlab = 'noise',
            #    ylab = '',
            #    type = c("o"),
            #    key = list(text=list(sort(ordered.regs.c5), col = mycolors),
            #            lines = list(lty = rep(1, length(ordered.regs.c5)), col = mycolors)),
            #    ylim = ylim))

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
                main = toupper(paste('mean edmonds vs scite sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c5), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c5)), col = mycolors)),
                ylim = ylim))
            #ggplot(res, aes(x=x, y=mean, colour=algorithm))
            #stop()
            dev.off()


            # edmonds vs scite 5 median

            res = results.c5[which(results.c5$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'edmods_scite_median', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'median']) - 0.05, max(res[,'median']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c5), 'Set1')[1:length(ordered.regs.c5)]
            print(xyplot(median~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('median edmonds vs scite sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c5), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c5)), col = mycolors)),
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
