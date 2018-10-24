plot_Spectrum2 <- function(spectrum,
                           reporters = NULL,
                           full = FALSE,
                           centroided. = centroided(spectrum),
                           plot = TRUE) {
    alpha <- NULL # to satisfy codetools 'no visible binding...'
    if (is.null(reporters)) {
        full <- TRUE
    } else {
        if (class(reporters) != "ReporterIons")
            stop("Reporters must be of class 'ReporterIons'.")
    }
    xmin <- xmax <- ymin <- ymax <- fill <- NULL # to satisfy codetools
    ## plot_Spectrum2: no visible binding for global variable 'xmin'
    ## ...
    mtc <- mz(spectrum)
    i <- intensity(spectrum)
    ## preparing full spectrum plot p
    dfr <- data.frame(i = i, mtc = mtc)
    if (nrow(dfr) == 0)
        stop("No data to be plotted in full scan")
    mainvp <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)
    title <- ggtitle(paste("Precursor M/Z",
                           round(precursorMz(spectrum), 2)))
    if (centroided.) {
        p <- ggplot(dfr, aes(x = mtc, xend = mtc, y = 0, yend = i)) +
            geom_segment() +
            theme(legend.position = "none") +
            labs(x = "M/Z", y = "Intensity")
    } else {
        p <- ggplot(dfr, aes(x = mtc, y = i)) + theme(legend.position = "none") +
            labs(x = "M/Z", y = "Intensity") +
            geom_line() ## + geom_point(alpha=I(1/10))
    }
    ## data reporters plot reps
    if (!is.null(reporters)) {
        width <- reporters@width
        rlim1 <- min(reporters@mz) - width
        rlim2 <- max(reporters@mz) + width
        dfr2 <- subset(dfr, mtc >= rlim1 & mtc <= rlim2)
        if ( nrow(dfr2) == 0 ) {
            warning("No reporter peaks to be plotted.")
            reporters <- NULL
        } else {
            coord <- coord_cartesian(xlim = c(rlim1, rlim2))
            subvp <- viewport(width = 2/3,
                              height = 1/3,
                              x = .95,
                              y = .90, ## was .92
                              just = c("right","top"))
            rectdfr <- data.frame(mtc = mean(dfr$mtc),
                                  i = 0,
                                  xmin = reporters@mz-reporters@width,
                                  xmax = reporters@mz+reporters@width,
                                  ymin = 0,
                                  ymax = max(dfr$i),
                                  fill = reporters@col,
                                  alpha = 0.2) ## was 1/3
            rect <- geom_rect(data = rectdfr,
                              aes(xmin = xmin,
                                  xmax = xmax,
                                  ymin = ymin,
                                  ymax = ymax,
                                  fill = fill,
                                  alpha = alpha))
            if (centroided.) {
                p2 <- ggplot(dfr2, aes(x = mtc, xend = mtc, y = 0, yend = i)) +
                    geom_segment() + rect
            } else {
                p2 <- ggplot(dfr2, aes(x = mtc, y = i)) + geom_line() + rect
            }
            reps <- p2 + coord +
                theme_gray(5) +
                labs(x = NULL, y = NULL) +
                theme(plot.margin = unit(c(1,1,0,0), "lines")) +
                scale_x_continuous(breaks=seq(rlim1, rlim2, (rlim2-rlim1)/10)) +
                theme(legend.position = "none")
        }
    }
    ## plotting
    if (full) {
        if (plot) {
            print(p + title, vp = mainvp)
            if ( !is.null(reporters) )
                print(reps, vp = subvp)
        }
        invisible(p + title)
    } else {
        if ( plot & !is.null(reporters) ) {
            print(reps + title, vp = mainvp)
            invisible(reps + title)
        }
    }
}
