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
                             verbose = isMSnbaseVerbose()) {
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
    delta <- data.frame(value = delta)
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

##' A function to convert the identification data contained in an
##' \code{mzRident} object to a \code{data.frame}. Each row represents
##' a scan, which can however be repeated several times if the PSM
##' matches multiple proteins and/or contains two or more
##' modifications. To reduce the \code{data.frame} so that rows/scans
##' are unique and use semicolon-separated values to combine
##' information pertaining a scan, use
##' \code{\link[=reduce,data.frame-method]{reduce}}.
##'
##' See also the \emph{Tandem MS identification data} section in the
##' \emph{MSnbase-demo} vignette.
##'
##' @title Coerce identification data to a \code{data.frame}
##' @param from An object of class \code{mzRident} defined in the
##'     \code{mzR} package.
##' @return A \code{data.frame}
##' @author Laurent Gatto
##' @name as
##' @rdname mzRident2dfr
##' @aliases as.data.frame.mzRident
##' @examples
##' ## find path to a mzIdentML file
##' identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
##'                  full.name = TRUE, pattern = "dummyiTRAQ.mzid")
##' library("mzR")
##' x <- openIDfile(identFile)
##' x
##' as(x, "data.frame")
setAs("mzRident", "data.frame",
      function(from) {
          ## peptide spectrum matching
          iddf <- factorsAsStrings(psms(from))
          ## add file raw and mzid provenances
          src <- basename(sourceInfo(from))
          if (length(src) > 1) ## see issue #261
              src <- paste(src, collapse = ";")
          iddf$spectrumFile <- src
          iddf$idFile <- basename(fileName(from))
          ## add scores
          scores <- factorsAsStrings(score(from))
          if (nrow(scores)) { ## see issue #261
              stopifnot(identical(iddf[, 1], scores[, 1]))
              iddf <- cbind(iddf, scores[, -1])
          }
          ## add modification
          mods <- factorsAsStrings(modifications(from))
          names(mods)[-1] <- makeCamelCase(names(mods), prefix = "mod")[-1]
          iddf <- merge(iddf, mods,
                        by.x = c("spectrumID", "sequence"),
                        by.y = c("spectrumID",  "modSequence"),
                        suffixes = c("", ".y"),
                        all = TRUE, sort = FALSE)
          iddf[, "spectrumID.y"] <- NULL
          ## add substitutions
          subs <- factorsAsStrings(substitutions(from))
          names(subs)[-1] <- makeCamelCase(names(subs), prefix = "sub")[-1]
          iddf <- merge(iddf, subs,
                        by.x = c("spectrumID" = "sequence"),
                        by.y = c("spectrumID" = "subSequence"),
                        suffixes = c("", ".y"),
                        all = TRUE, sort = FALSE)
          iddf[, "spectrumID.y"] <- NULL
          iddf
      })

as.data.frame.mzRident <-
    function(x, row.names = NULL, optional = FALSE, ...) as(x, "data.frame")
