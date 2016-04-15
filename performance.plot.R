library(lattice)
library(RColorBrewer)

performance.plot <- function(dataset,
    sample.type,
    branching) {

    if (! branching %in% c('low', 'medium', 'high')) {
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

                zlim = range(stats)
                zlim[1] = zlim[1] - 0.05
                zlim[2] = 1.0

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

    if (! branching %in% c('low', 'medium', 'high')) {
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

        algorithms = get(type, dataset)
        for (algorithm in names(algorithms)) {
            cat('\n  ', algorithm, ' ')
            regs = get(algorithm, algorithms)
            for (reg in names(regs)) {
                if ((algorithm == 'caprese' && reg == 'no.reg')
                    || (algorithm == 'prim' && reg == 'no.reg')
                    || (algorithm == 'edmonds' && reg == 'no.reg')
                    || (algorithm == 'chowliu' && reg == 'loglik')) {
                    #|| (algorithm == 'capri' && reg == 'bic')) {

                    cat('\n    ', reg)
                    stats = get(reg, regs)
                    for (i in 1:nrow(stats)) {
                        for(j in 1:ncol(stats)) {
                            results = rbind(results, c(as.numeric(i),
                                as.numeric(j),
                                as.numeric(stats[i,j]),
                                algorithm), stringsAsFactors = FALSE)
                        }
                    }
                    ordered.regs = c(ordered.regs, algorithm)
                    next
                }
            }
        }

        colnames(results) = c('x', 'y', 'z', 'algorithm')
        results$x = as.numeric(results$x)
        results$y = as.numeric(results$y)
        results$z = as.numeric(results$z)
        results$algorithm = as.factor(results$algorithm)

        pdf(file=paste0(branching, '/', sample.type, '/', type, '/', 'compare', '.pdf'), width=8.5, height=6.5)

        zlim = c(min(results[,'z']) - 0.05, max(results[,'z']) + 0.05)
        mycolors = brewer.pal(length(ordered.regs), 'Accent')
        mycolors.trans = add.alpha(mycolors, 0.7)

        print(wireframe(z~x*y,
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
        cat('\n')
    }
}


compare.performance.plot.2d <- function(dataset,
    sample.type,
    branching) {

    if (! branching %in% c('low', 'medium', 'high')) {
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
        ordered.regs.c1 = NULL
        ordered.regs.c2 = NULL
        ordered.regs.c3 = NULL
        ordered.regs.c4 = NULL

        algorithms = get(type, dataset)
        for (algorithm in names(algorithms)) {
            cat('\n  ', algorithm, ' ')
            regs = get(algorithm, algorithms)
            for (reg in names(regs)) {

                # cfg 1 - mix
                if ((algorithm == 'caprese' && reg == 'no.reg')
                    || (algorithm == 'prim' && reg == 'no.reg')
                    || (algorithm == 'edmonds' && reg == 'no.reg')
                    || (algorithm == 'chowliu' && reg == 'loglik')
                    || (algorithm == 'capri' && reg == 'bic')) {

                    cat('\n    ', reg)
                    stats = get(reg, regs)
                    for (i in 1:nrow(stats)) {
                        for(j in 1:ncol(stats)) {
                            if (sample[j] %in% select.sample) {
                                results.c1 = rbind(results.c1, c(i,
                                    stats[i,j],
                                    sample[j],
                                    algorithm), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c1 = c(ordered.regs.c1, algorithm)
                }

                # cfg 2 - all prim
                if (algorithm == 'prim') {

                    cat('\n    ', reg)
                    stats = get(reg, regs)
                    for (i in 1:nrow(stats)) {
                        for(j in 1:ncol(stats)) {
                            if (sample[j] %in% select.sample) {
                                results.c2 = rbind(results.c2, c(i,
                                    stats[i,j],
                                    sample[j],
                                    reg), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c2 = c(ordered.regs.c2, reg)
                }

                # cfg 3 - all edmonds
                if (algorithm == 'edmonds') {
                    cat('\n    ', reg)
                    stats = get(reg, regs)
                    for (i in 1:nrow(stats)) {
                        for(j in 1:ncol(stats)) {
                            if (sample[j] %in% select.sample) {
                                results.c3 = rbind(results.c3, c(i,
                                    stats[i,j],
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
                    stats = get(reg, regs)
                    for (i in 1:nrow(stats)) {
                        for(j in 1:ncol(stats)) {
                            if (sample[j] %in% select.sample) {
                                results.c4 = rbind(results.c4, c(i,
                                    stats[i,j],
                                    sample[j],
                                    reg), stringsAsFactors = FALSE)
                            }
                        }
                    }
                    ordered.regs.c4 = c(ordered.regs.c4, reg)
                }

            }
        }

        colnames(results.c1) = c('x', 'y', 'sample', 'algorithm')
        colnames(results.c2) = c('x', 'y', 'sample', 'algorithm')
        colnames(results.c3) = c('x', 'y', 'sample', 'algorithm')
        colnames(results.c4) = c('x', 'y', 'sample', 'algorithm')
        results.c1$x = as.numeric(results.c1$x)
        results.c1$y = as.numeric(results.c1$y)
        results.c1$sample = as.numeric(results.c1$sample)
        results.c1$algorithm = as.factor(results.c1$algorithm)
        results.c2$x = as.numeric(results.c2$x)
        results.c2$y = as.numeric(results.c2$y)
        results.c2$sample = as.numeric(results.c2$sample)
        results.c2$algorithm = as.factor(results.c2$algorithm)
        results.c3$x = as.numeric(results.c3$x)
        results.c3$y = as.numeric(results.c3$y)
        results.c3$sample = as.numeric(results.c3$sample)
        results.c3$algorithm = as.factor(results.c3$algorithm)
        results.c4$x = as.numeric(results.c4$x)
        results.c4$y = as.numeric(results.c4$y)
        results.c4$sample = as.numeric(results.c4$sample)
        results.c4$algorithm = as.factor(results.c4$algorithm)


        for (s in select.sample) {

            # prim res 1

            res = results.c1[which(results.c1$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'mix', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'y']) - 0.05, max(res[,'y']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c1), 'Dark2')
            print(xyplot(y~x,
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

            res = results.c2[which(results.c2$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'prim', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'y']) - 0.05, max(res[,'y']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c2), 'Set1')
            print(xyplot(y~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('prim sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c2), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c2)), col = mycolors)),
                ylim = ylim))

            dev.off()


            # edmonds res 3

            res = results.c3[which(results.c3$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'edmonds', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'y']) - 0.05, max(res[,'y']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c3), 'Set1')
            print(xyplot(y~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('edmonds sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c3), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c3)), col = mycolors)),
                ylim = ylim))

            dev.off()

            # chowliu res 4

            res = results.c4[which(results.c4$sample == s), ]
            pdf(file = paste0(branching, '/', sample.type, '/', type, '/', s, '_', 'chowliu', '.pdf'), width=8.5, height=6.5)
            ylim = c(min(res[,'y']) - 0.05, max(res[,'y']) + 0.05)
            mycolors = brewer.pal(length(ordered.regs.c4), 'Set1')
            print(xyplot(y~x,
                grid = TRUE,
                data = res,
                groups = algorithm,
                col = mycolors,
                main = toupper(paste('chowliu sample:', s, branching, type)),
                xlab = 'noise',
                ylab = '',
                type = c("o"),
                key = list(text=list(sort(ordered.regs.c4), col = mycolors),
                        lines = list(lty = rep(1, length(ordered.regs.c4)), col = mycolors)),
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
