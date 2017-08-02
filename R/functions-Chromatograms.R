## functions for Chromatograms class

.validChromatograms <- function(x) {
    msg <- character()
    ## All elements have to be of type Chromatogram
    if (length(x)) {
        res <- vapply(x, FUN = is, FUN.VALUE = logical(1L), "Chromatogram")
        if (!all(res))
            msg <- c(msg, paste0("All elements have to be of type ",
                                 "'Chromatogram'."))
        ## Shall we also ensure that fromFile in each column is the same?
    }
    if (nrow(x@phenoData) != ncol(x))
        msg <- c(msg, paste0("nrow of phenoData has to match ncol ",
                             "of the Chromatograms object"))
    ## Check colnames .Data with rownames phenoData.
    if (!is.null(colnames(x))) {
        if (any(colnames(x) != rownames(x@phenoData)))
            msg <- c(msg, paste0("colnames of object has to match rownames of",
                                 " phenoData"))
    } else {
        if (any(rownames(x@phenoData) != as.character(1:ncol(x))))
            msg <- c(msg, paste0("rownames of phenoData does not match ",
                                 "colnames of object"))
    }
    if (length(msg))
        msg
    else TRUE
}

#' @description \code{Chromatograms}: create an instance of class
#'     \code{Chromatograms}.
#'
#' @param data A \code{list} of \code{\link{Chromatogram}} objects.
#'
#' @param phenoData either a \code{data.frame}, \code{AnnotatedDataFrame} or
#'     \code{NAnnotatedDataFrame} describing the phenotypical information of the
#'     samples.
#'
#' @param ... Additional parameters to be passed to the
#'     \code{\link[base]{matrix}} constructor, such as \code{nrow}, \code{ncol}
#'     and \code{byrow}.
#' 
#' @rdname Chromatograms-class
Chromatograms <- function(data, phenoData, ...) {
    if (missing(data))
        return(new("Chromatograms"))
    datmat <- matrix(data, ...)
    if (missing(phenoData))
        phenoData <- annotatedDataFrameFrom(datmat, byrow = FALSE)
    if (ncol(datmat) != nrow(phenoData))
        stop("Dimensions of the data matrix and the  phenoData do not match")
    ## Convert phenoData...
    if (is(phenoData, "data.frame"))
        phenoData <- AnnotatedDataFrame(phenoData)
    if (is(phenoData, "AnnotatedDataFrame"))
        phenoData <- as(phenoData, "NAnnotatedDataFrame")
    ## Set colnames if we have some
    if (any(rownames(phenoData) != as.character(1:nrow(phenoData))))
        colnames(datmat) <- rownames(phenoData)
    else
        rownames(phenoData) <- colnames(datmat)
    res <- new("Chromatograms", .Data = datmat, phenoData = phenoData)
    ## res@.Data <- datmat
    ## res@phenoData <- phenoData
    ## res <- as(matrix(data, ...), "Chromatograms")
    if (validObject(res))
        res
}


#' @description Plot the data from a list of Chromatogram objects (all
#'     representing the same MS data slice across multiple files) into the
#'     same plot.
#'
#' @note We are using the matplot here, since that is much faster than lapply
#'     on the individual chromatogram objects.
#'
#' @author Johannes Rainer
#'
#' @noRd
.plotChromatogramList <- function(x, col = "#00000060", lty = 1, type = "l",
                                  xlab = "retention time", ylab = "intensity",
                                  main = NULL, ...) {
    if (!is.list(x) & !all(vapply(x, FUN = is, FUN.VALUE = logical(1L),
                                  "Chromatogram")))
        stop("'x' has to be a list of Chromatogram objects")
    ## Check col, lty and type parameters
    if (length(col) != length(x))
        col <- rep(col[1], length(x))
    if (length(lty) != length(x))
        lty <- rep(lty[1], length(x))
    if (length(type) != length(x))
        type <- rep(type, length(x))
    if (is.null(main)) {
        suppressWarnings(
            mzr <- range(lapply(x, mz), na.rm = TRUE, finite = TRUE)
        )
        main <- paste0(format(mzr, digits = 7), collapse = " - ")
    }
    ## Number of measurements we've got per chromatogram. This can be different
    ## between samples, from none (if not a single measurement in the rt/mz)
    ## to the number of data points that were actually measured.
    lens <- lengths(x)
    maxLens <- max(lens)
    
    ints <- rts <- matrix(NA_real_, nrow = maxLens, ncol = length(x))
    for (i in seq(along = x)) {
        if (lens[i]) {
            rows <- seq_len(lens[i])
            rts[rows, i] <- rtime(x[[i]])
            ints[rows, i] <- intensity(x[[i]])
        }
    }
    ## Identify columns/samples that have only NAs in the intensity matrix.
    ## Such columns represent samples for which no valid intensity was measured
    ## in the respective mz slice (these would still have valid retention time
    ## values), or samples that don't have a single scan in the respective rt
    ## range.
    keep <- colSums(!is.na(rts)) > 0

    ## Finally plot the data.
    if (any(keep)) {
        matplot(x = rts[, keep, drop = FALSE],
                y = ints[, keep, drop = FALSE], type = type[keep],
                lty = lty[keep], col = col[keep], xlab = xlab,
                ylab = ylab, main = main, ...)
    } else {
        warning("No chromatographic data to plot")
        plot(x = 3, y = 3, pch = NA, xlab = xlab, ylab = ylab, main = main, ...)
    }
}


