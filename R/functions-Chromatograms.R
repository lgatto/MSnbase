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
    if (any(colnames(x) != rownames(x@phenoData)))
        msg <- c(msg, paste0("colnames of object has to match rownames of",
                             " object's phenoData"))
    if (nrow(x@featureData) != nrow(x))
        msg <- c(msg, paste0("nrow of featureData has to match nrow ",
                             "of the Chromatograms object"))
    if (any(rownames(x) != rownames(x@featureData)))
        msg <- c(msg, paste0("rownames of object has to match rownames of",
                             " object's featureData"))
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
#' @param featureData either a \code{data.frame} or \code{AnnotatedDataFrame}
#'     with additional information for each row of chromatograms.
#' 
#' @param ... Additional parameters to be passed to the
#'     \code{\link[base]{matrix}} constructor, such as \code{nrow}, \code{ncol}
#'     and \code{byrow}.
#' 
#' @rdname Chromatograms-class
Chromatograms <- function(data, phenoData, featureData, ...) {
    if (missing(data))
        return(new("Chromatograms"))
    datmat <- matrix(data, ...)
    if (missing(phenoData))
        phenoData <- annotatedDataFrameFrom(datmat, byrow = FALSE)
    if (ncol(datmat) != nrow(phenoData))
        stop("phenoData has to have the same number of rows as the data ",
             "matrix has columns")
    ## If colnames of datmat are NULL, use the rownames of phenoData
    if (is.null(colnames(datmat)))
        colnames(datmat) <- rownames(phenoData)
    ## Convert phenoData...
    if (is(phenoData, "data.frame"))
        phenoData <- AnnotatedDataFrame(phenoData)
    if (is(phenoData, "AnnotatedDataFrame"))
        phenoData <- as(phenoData, "NAnnotatedDataFrame")
    if (missing(featureData))
        featureData <- annotatedDataFrameFrom(datmat, byrow = TRUE)
    if (nrow(datmat) != nrow(featureData))
        stop("featureData has to have the same number of rows as the data ",
             "matrix has rows")
    if (is.null(rownames(datmat)))
        rownames(datmat) <- rownames(featureData)
    if (is(featureData, "data.frame"))
        featureData <- AnnotatedDataFrame(featureData)
    res <- new("Chromatograms", .Data = datmat, phenoData = phenoData,
               featureData = featureData)
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
    keep <- colSums(!is.na(ints)) > 0

    ## Finally plot the data.
    if (any(keep)) {
        matplot(x = rts[, keep, drop = FALSE],
                y = ints[, keep, drop = FALSE], type = type[keep],
                lty = lty[keep], col = col[keep], xlab = xlab,
                ylab = ylab, main = main, ...)
    } else {
        warning("Chromatograms empty")
        plot(3, 3, pch = NA, xlab = xlab, ylab = ylab, main = main)
        text(3, 3, labels = "Empty Chromatograms", col = "red")
    }
}

#' Helper function to extract mz, precursorMz or productMz from a Chromatograms
#' object
#'
#' @author Johannes Rainer
#'
#' @noRd
.mz_chromatograms <- function(x, mz = c("mz", "precursorMz", "productMz")) {
    mz <- match.arg(mz)
    if (!nrow(x))
        return(matrix(nrow = 0, ncol = 2, dimnames = list(character(),
                                                          c("mzmin", "mzmax"))))
    ## If we've got the values in the featureData, use these.
    if (mz %in% c("precursorMz", "productMz"))
        vl <- rep(sub("Mz", "IsolationWindowTargetMZ", mz), 2)
    else
        vl <- c("mzmin", "mzmax")
    if (all(vl %in% fvarLabels(x))) {
        ## Want to return a matrix, not a data.frame
        cbind(mzmin = fData(x)[, vl[1]], mzmax = fData(x)[, vl[2]])
    } else {
        ## Get the xxx mz from the Chromatogram objects. Throw an error if
        ## the values in one row are not identical
        mzr <- matrix(nrow = nrow(x), ncol = 2,
                      dimnames = list(NULL, c("mzmin", "mzmax")))
        for (i in seq_len(nrow(mzr))) {
            rngs <- unique(do.call(
                rbind, lapply(x@.Data[i, ], getMethod(mz, "Chromatogram"))))
            if (nrow(rngs) != 1)
                stop("Chromatograms in row ", i, " have different ", mz)
            mzr[i, ] <- rngs
        }
        mzr
    }
}

#' Simple binning function for Chromatograms object. Defines common breaks for
#' `Chromatogram` objects in each row.
#'
#' @author Johannes Rainer
#' 
#' @noRd
.bin_Chromatograms <- function(object, binSize = 0.5, breaks = numeric(),
                               fun = max) {
    for (i in seq_len(nrow(object))) {
        if (!length(breaks)) {
            rt_rng <- range(lapply(object[i, ], function(z) range(rtime(z))))
            brks <- .fix_breaks(seq(floor(rt_rng[1]), ceiling(rt_rng[2]),
                                    by = binSize), rt_rng)
        } else brks <- breaks
        object[i, ] <- lapply(object[i, ], .bin_Chromatogram, binSize = binSize,
                              breaks = brks, fun = fun)
    }
    if (validObject(object))
        object
}

