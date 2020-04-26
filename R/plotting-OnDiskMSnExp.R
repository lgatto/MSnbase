## Test whether we could speed up things for OnDiskMSnExp objects.
plot_OnDiskMSnExp <- function(object,
                              reporters,
                              full = FALSE,
                              centroided.,
                              plot = TRUE) {
    i <- NULL # to satisfy codetools
    ## plot_MSnExp: no visible binding for global variable i
    if (missing(centroided.))
        centroided. <- any(centroided(object))
    mtc <- unlist(mz(object))
    spectraList <- spectra(object)
    ints <- unlist(sapply(spectraList, function(x) x@intensity))
    mzs <- unlist(sapply(spectraList, function(x) x@mz))
    l <- unlist(sapply(spectraList, function(x) length(x@mz)))
    n <- rep(1:length(l), l)
    dfr <- data.frame(i = ints, mz = mzs, n = n)
    colnames(dfr) <- c("i", "mz", "n")
    if (all(msLevel(object) > 1)) {
        pmz <- paste(unique(unlist(sapply(spectraList,
                                          function(x) round(precursorMz(x), 2)))),
                     collapse = ",")
        title <- ggtitle(paste("Precursor M/Z", pmz))
    } else {
        rtm <- paste(formatRt(range(rtime(object))), collapse = " - ")
        title <- ggtitle(paste("Retention time", rtm))
        full <- TRUE
    }
    if (centroided.) {
        p <- ggplot(dfr, aes(x = mz, xend = mz, y = 0, yend = i)) +
            geom_segment()
    } else {
        p <- ggplot(data = dfr, aes(x = mz, y = i)) + geom_line()
    }
    p <- p +  facet_grid(n~.) +
        labs(x = "M/Z", y = "Intensity") +
            title
    if (!full) {
        if (class(reporters) != "ReporterIons")
            stop("Reporters must be of class \"ReporterIons\".")
        width <- reporters@width
        rlim1 <- min(reporters@mz) - width
        rlim2 <- max(reporters@mz) + width
        reps <- coord_cartesian(xlim = c(rlim1, rlim2))
        breaks <- scale_x_continuous(breaks = seq(rlim1, rlim2, (rlim2-rlim1)/10))
        p <- p + reps + breaks
    }
    if (plot)
        print(p)
    invisible(p)
}


