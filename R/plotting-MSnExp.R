plot_MSnExp <- function(object,
                        reporters,
                        full = FALSE,
                        centroided.,
                        plot = TRUE) {
    i <- NULL # to satisfy codetools
    ## plot_MSnExp: no visible binding for global variable 'i'
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
    p <- p +
        facet_grid(n ~ ., scales = "free_y") +
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

plotMzDelta_MSnExp <- function(object,            ## MSnExp object
                               reporters = NULL,  ## reporters to be removed
                               percentage = 0.1,  ## percentage of peaks to consider
                               precMz = NULL,     ## precursors to be removed
                               precMzWidth = 2,   ## precrsor m/z with
                               bw = 1,            ## histogram bandwidth
                               xlim = c(40,200),  ## delta m/z range
                               withLabels = TRUE, ## add amino acide labels
                               size = 2.5,        ## labels size
                               plot = TRUE,       ## plot figure
                               verbose = isMSnbaseVerbose()) {
    ## Contributed by Guangchuang Yu for the plotMzDelta QC
    ## Modified aa labelling
    ResidueMass <- ..density.. <- NULL ## to accomodate codetools
    value <- AA <- NULL
    delta <- vector("list", length(object))
    if (verbose) {
        pb <- txtProgressBar(min = 1, max = length(object), style = 3)
        k <- 1
    }
    if (!is.null(reporters)) {
        if (verbose)
            message("Removing reporter ion peaks...\n")
        object <- removeReporters(object, reporters, verbose = FALSE)
    }
    ## TODO -- check if there is a precursor mz to remove,
    ## i.e precursorMz != 0
    for (j in 1:length(object)) {
        if (verbose) {
            setTxtProgressBar(pb, k)
            k <- k + 1
        }
        sp <- object[[j]]
        ## TODO - better than setting precMzWidth statically
        ## would be to get the peaks based on it m/z value
        ## and then find it's upper/lower m/z limits to set to 0
        sp <- utils.removePrecMz_Spectrum(sp)
        ds <- utils.getMzDelta(sp, percentage)
        delta[[j]] <- ds[ds > xlim[1] & ds < xlim[2]]
    }
    if (verbose) {
        close(pb)
        message(" Plotting...\n")
    }
    delta <- data.frame(value = unlist(delta))
    p <- ggplot(delta, aes(x = value)) +
        geom_histogram(aes(y = ..density..), stat = "bin", binwidth = bw) +
            scale_x_continuous(limits = xlim) +
                    xlab("m/z delta") + ylab("Density") +
                        ggtitle("Histogram of Mass Delta Distribution")
    if (withLabels) {
        y_offset <- x_offset <- rep(0.5, 21)
        names(y_offset) <- names(x_offset) <- PSMatch::getAminoAcids()$AA
        x_offset[c("I", "L", "K", "Q")] <- 1
        y_offset[c("V", "C")] <- 1
        y_offset[c("P", "T")] <- 0
        y_offset[c("N", "E")] <- 1
        y_offset[c("K", "Q", "I", "L")] <- 0
        y_offset[c("D", "M")] <- 0
        aa <- cbind(PSMatch::getAminoAcids(), x_offset, y_offset)
        ## removing Isoleucine, as it has the same residue mass
        ## as leucine, and updating leucine's label to I/L
        aa$AA <- as.character(aa$AA)
        aa[aa$AA == "I", "ResidueMass"] <- NA
        aa[aa$AA == "L", "AA"] <- "I/L"
        ## Removing Q as it is too close to K to show
        ## up correctly and updating K label to K/Q
        aa[aa$AA == "Q", "ResidueMass"] <- NA
        aa[aa$AA == "K", "AA"] <- "K/Q"
        p <- p +
            geom_vline(data = aa,
                       aes(xintercept = ResidueMass,
                           colour = AA),
                       alpha = I(1/2)) +
                           geom_text(data = aa,
                                     aes(x = ResidueMass,
                                         y = -0.001, label = AA,
                                         vjust = y_offset,
                                         hjust = x_offset),
                                     size = size) +
                                         theme(legend.position = "none")
    }
    if (plot)
        print(p)
    invisible(p)
}
