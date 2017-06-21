## functions for Chromatograms class

.validChromatograms <- function(x) {
    msg <- character()
    ## All elements have to be of type Chromatogram
    if (length(x)) {
        res <- sapply(x, FUN = function(z) is(z, "Chromatogram"))
        if (!all(res))
            msg <- c(msg, paste0("All elements have to be of type ",
                                 "'Chromatogram'."))
        ## Shall we also ensure that fromFile in each column is the same?
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
#' @param ... Additional parameters to be passed to the \code{\link{matrix}}
#'     constructor, such as \code{nrow}, \code{ncol} and \code{byrow}.
#' 
#' @rdname Chromatograms-class
Chromatograms <- function(data, ...) {
    if (missing(data))
        return(new("Chromatograms"))
    res <- as(matrix(data, ...), "Chromatograms")
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
    if (!is.list(x) & !all(sapply(x, function(z) is(z, "Chromatogram"))))
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
    lens <- unique(lengths(x))
    max_len <- max(lens)
    max_len_vec <- rep_len(NA, max_len)
    ## Generate the matrix of rt values, columns are samples, rows retention
    ## time values. Fill each column with NAs up to the maximum number of values
    ## we've got in a sample/file.
    rts <- do.call(cbind, lapply(x, function(z) {
        cur_len <- length(z)
        if (cur_len == 0)
            max_len_vec
        else {
            max_len_vec[seq_len(cur_len)] <- rtime(z)
            max_len_vec
        }
    }))
    ## Same for the intensities.
    ints <- do.call(cbind, lapply(x, function(z) {
        cur_len <- length(z)
        if (length(z) == 0)
            max_len_vec
        else {
            max_len_vec[seq_len(cur_len)] <- intensity(z)
            max_len_vec
        }
    }))
    ## Identify columns that have only NAs in either intensity or rt - these
    ## will not be plotted.
    keepCol <- which(apply(ints, MARGIN = 2, function(z) any(!is.na(z))) |
                     apply(rts, MARGIN = 2, function(z) any(!is.na(z))))
    ## Finally plot the data.
    if (length(keepCol)) {
        matplot(x = rts[, keepCol, drop = FALSE],
                y = ints[, keepCol, drop = FALSE], type = type[keepCol],
                lty = lty[keepCol], col = col[keepCol], xlab = xlab,
                ylab = ylab, main = main, ...)
    } else
        plot(x = 3, y = 3, pch = NA, xlab = xlab, ylab = ylab, main = main, ...)
}
