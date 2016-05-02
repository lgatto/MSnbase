peaksAsLists <-  function(object,
                          i,
                          what = c("both", "mz", "int"),
                          simplify = FALSE) {
    what <- match.arg(what)
    if (missing(i)) {
        pl <- peaks(object)
    } else {
        pl <- peaks(object, i)
        if (is.matrix(pl)) ## i was of length 1
            pl <- list(pl)
    }
    switch(what,
           mz = lapply(pl, function(x) list(mz = x[, 1])),
           int = lapply(pl, function(x) list(int = x[, 2])),
           both = lapply(pl, function(x) list(mz = x[, 1],
                                              int = x[, 2])))
}


list2Spectrum2 <- function(x, ...)
    new("Spectrum2",
        intensity = x$int,
        mz = x$mz,
        ...)

plotMzDelta_list <- function(object,            ## peakLists
                             reporters = NULL,  ## reporters to be removed
                             percentage = 0.1,  ## percentage of peaks to consider
                             precMz,            ## precursors to be removed
                             precMzWidth = 2,   ## precrsor m/z with
                             bw = 1,            ## histogram bandwidth
                             xlim = c(40, 200), ## delta m/z range
                             withLabels = TRUE, ## add amino acide labels
                             size = 2.5,        ## labels size
                             plot = TRUE,       ## plot figure
                             verbose = TRUE) {
    if (missing(precMz))
        stop("Precursor M/Z is only supported for MSnExp instances.")
    ResidueMass <- ..density.. <- NULL ## to accomodate codetools
    value <- AA <- NULL
    delta <- vector("list", length(object))
    if (verbose) {
        pb <- txtProgressBar(min = 1, max = length(object), style = 3)
        k <- 1
    }
    if (!is.null(reporters)) {
        warning("Currenlty only available for MSnExp instances.")
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
        sp <- utils.removePrecMz_list(sp, precMz[j], precMzWidth)
        ds <- utils.getMzDelta_list(sp, percentage)
        delta[[j]] <- ds[ds > xlim[1] & ds < xlim[2]]
    }   
    if (verbose) {
        close(pb)
        message(" Plotting...\n")
    }
    delta <- unlist(delta)
    ## could round deltas to speed up?
    delta <- melt(delta)
    p <- ggplot(delta, aes(x = value)) + 
        geom_histogram(aes(y = ..density..), stat = "bin", binwidth = bw) +
            scale_x_continuous(limits = xlim) +
                    xlab("m/z delta") + ylab("Density") +
                        ggtitle("Histogram of Mass Delta Distribution")
    if (withLabels) {
        y_offset <- x_offset <- rep(0.5, 21)
        names(y_offset) <- names(x_offset) <- .get.amino.acids()$AA
        x_offset[c("I", "L", "K", "Q")] <- 1
        y_offset[c("V", "C")] <- 1
        y_offset[c("P", "T")] <- 0
        y_offset[c("N", "E")] <- 1
        y_offset[c("K", "Q", "I", "L")] <- 0
        y_offset[c("D", "M")] <- 0
        aa <- cbind(.get.amino.acids(), x_offset, y_offset)
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


.chromatogram <- function(hd,
                          y = c("tic", "bpi"),
                          f,
                          legend = TRUE,
                          plot = TRUE,
                          ...) {
    y <- match.arg(y)
    ylab <- switch(y,
                tic = "Total ion current",
                bpi = "Base peak intensity")
    yy <- switch(y,
                tic = hd$totIonCurrent,
                bpi = hd$basePeakIntensity)
    yymax <- max(yy)
    yy <- yy / yymax * 100
    xx <- hd$retentionTime
    if (plot) {
        plot(yy ~ xx, type = "l",
             xlab = "Time (sec)", ylab = ylab,
             ...)
        abline(h = 0)
        if (legend) {
            leg <- sprintf("%s: %.3g", toupper(y), yymax)
            if (!missing(f))
                leg <- c(f, leg)
            legend("topleft",
                   leg,
                   col = NA,
                   cex = .7,
                   bty = "n")
        }
    }
    dd <- data.frame(xx, yy)
    colnames(dd) <- c("rt", y)
    invisible(dd)
}

xicplot <- function(dd, mz, width, rtlim,
                    npeaks, legend, points, ...) {
    ptcex <- .5
    plot(int ~ rt, data = dd,
         ylab = "Counts", xlab = "Retention time [s]",
         type = "l", xlim = rtlim,
         ...)
    abline(h = 0, col = "grey")
    grid()
    if (legend)
        legend("topleft",
               c(paste0("Precursor: ", mz),
                 paste0("XIC: ", mz-width, " - ", mz+width)),
               bty = "n", cex = .75)
    if (npeaks > 0) {
        dd2 <- dd[dd$rt >= min(rtlim) & dd$rt <= max(rtlim), ]
        dd2$int[dd2$int < max(dd2$int)/100] <- 0
        kx <- ky <- kz<- rep(NA, npeaks)
        for (.k in 1:length(kx)) {
            i <- order(dd2$int, decreasing = TRUE)[1]
            kx[.k] <- dd2$rt[i]
            ky[.k] <- dd2$int[i]
            kz[.k] <- dd2$mz[i]
            dd2$int[i] <- 0
            dd0 <- which(dd2$int == 0)
            ii <- which(dd0 == i)
            if (length(ii) < 3) break
            dd2$int[dd0[ii - 1]:dd0[ii + 1]] <- 0
        }
        ksel <- !is.na(kx)
        text(kx[ksel],
             ky[ksel],
             sprintf("%.4f", kz[ksel]),
             pos = 3,
             cex = .75)
    }
    if (points) {
        ## relevant MS2 spectra are coloured in xic_1
        points(int ~ rt, data = dd,
               col = "#00000060",
               pch = 19, cex = ptcex)
        ## highlight annotated peaks
        points(kx, ky, cex = ptcex, pch = 19,
               col = "#FFA404FF")
    }
}


## want a vectorised version
xic_1 <- function(object, ##
                  mz,
                  width = 0.5,
                  rtlim,
                  npeaks = 3,
                  charge,
                  clean = TRUE,
                  legend=TRUE,
                  plot = TRUE,
                  points = TRUE,
                  hd,
                  ...) {
    if (!missing(charge))
        mz <- mz/as.integer(charge)
    if (missing(hd))
        hd <- header(object)
    ms1 <- which(hd$msLevel == 1)
    pl <- peaks(object, ms1)
    hd1 <- hd[ms1, ]
    res <- lapply(seq_len(length(pl)),
                  function(i) {
                      ## matching peaks
                      sel <- pl[[i]][, 1] > mz - width &
                          pl[[i]][, 1] < mz + width
                      ans <- NA
                      if (any(sel)) {
                          j <- which.max(pl[[i]][sel, 2])
                          ## index, intensity, mz
                          ans <- c(i, pl[[i]][sel, 2][j],
                                   pl[[i]][sel, 1][j])
                      }
                      ans
                  })
    if (length(res2 <- res[!is.na(res)]) == 0)
        stop("No matching peaks found.")
    if (length(res2) < 15)
        warning("Only ", length(res2), " matching spectra found.")
    dd <- data.frame(int = sapply(res2, "[", 2),
                     rt = hd1$retentionTime[sapply(res2, "[", 1)],
                     mz = sapply(res2, "[", 3))
    if (clean)
        dd <- dd[utils.clean(dd$int, all = FALSE), ]
    if (plot) {
        if (missing(rtlim))
            rtlim <- range(dd$rt)
        xicplot(dd, mz, width, rtlim, npeaks, legend, points, ...)
        if (points) {
            hd2 <- hd[-ms1, ]
            sel <- hd2$precursorMZ > mz - width &
                hd2$precursorMZ < mz + width
            if (any(sel)) {
                pi <- hd2[sel, "precursorScanNum"]
                pj <- match(pi, hd1$acquisitionNum)
                ## pl2 <- peaks(object, pj) ## relevant MS2 spectra
                .dd2 <- data.frame(int = sapply(res[pj], "[", 2),
                                   rt = hd1$retentionTime[sapply(res[pj], "[", 1)],
                                   mz = sapply(res[pj], "[", 3))
                points(.dd2$rt, .dd2$int, col = "#FF0000",
                       pch = 19, cex = .6)
            }
        }
    }
    invisible(dd)
}
