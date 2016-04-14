library(lattice)

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