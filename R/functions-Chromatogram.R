.SUPPORTED_AGG_FUN_CHROM <- c("sum", "max", "min", "mean")
names(.SUPPORTED_AGG_FUN_CHROM) <-
    c("Total ion chromatogram (TIC).", "Base peak chromatogram (BPC).",
      "Intensity representing the minimum intensity across the mz range.",
      "Intensity representing the mean intensity across the mz range.")

#' @title Validation function for Chromatogram objects
#'
#' @description This function can be used instead of the \code{validObject} to
#'     check if the chromatogram is valid, without having to call the validity
#'     method on all super classes.
#'
#' @param object A \code{Chromatogram} object.
#'
#' @return \code{TRUE} if the \code{object} is valid and the error messages
#'     otherwise (i.e. a \code{character}).
#'
#' @author Johannes Rainer
#'
#' @noRd
.validChromatogram <- function(object) {
    msg <- character()
    if (length(object@rtime) != length(object@intensity))
        msg <- c(msg, "Length of 'rt' and 'intensity' have to match!")
    if (is.unsorted(object@rtime))
        msg <- c(msg, paste0("'rtime' has to be increasingly ordered!"))
    if (length(object@mz) > 0 & length(object@mz) != 2)
        msg <- c(msg, paste0("'mz' is supposed to contain the ",
                             "minimum and maximum mz values for the ",
                             "chromatogram."))
    if (!all(is.na(object@mz)))
        if (is.unsorted(object@mz))
            msg <- c(msg, "'mz' has to be increasingly ordered!")
    if (length(object@filterMz) > 0 & length(object@filterMz) != 2)
        msg <- c(msg, paste0("'filterMz' is supposed to contain the ",
                             "minimum and maximum mz values of the filter",
                             " used to create the chromatogram."))
    if (length(object@precursorMz) > 0 & length(object@precursorMz) != 2)
        msg <- c(msg, paste0("'precursorMz' is supposed to be a numeric of",
                             " length 2."))
    if (length(object@productMz) > 0 & length(object@productMz) != 2)
        msg <- c(msg, paste0("'productMz' is supposed to be a numeric of",
                             " length 2."))
    if (length(object@fromFile) > 1 | any(object@fromFile < 0))
        msg <- c(msg, paste0("'fromFile' is supposed to be a single ",
                             "positive integer!"))
    if (length(object@aggregationFun) > 1)
        msg <- c(msg, "Length of 'aggregationFun' has to be 1!")
    if (length(object@aggregationFun)) {
        if (!object@aggregationFun %in% .SUPPORTED_AGG_FUN_CHROM)
            msg <- c(msg, paste0("Invalid value for 'aggregationFun'! only ",
                                 paste0("'", .SUPPORTED_AGG_FUN_CHROM,"'",
                                        collapse = ","), " are allowed!"))
    }
    if (length(msg))
        msg
    else
        TRUE
}

#' @rdname Chromatogram-class
Chromatogram <- function(rtime = numeric(), intensity = numeric(),
                         mz = c(NA_real_, NA_real_),
                         filterMz = c(NA_real_, NA_real_),
                         precursorMz = c(NA_real_, NA_real_),
                         productMz = c(NA_real_, NA_real_),
                         fromFile = integer(),
                         aggregationFun = character(),
                         msLevel = 1L) {
    ## Check if we have to re-order the data (issue #145).
    if (is.unsorted(rtime)) {
        idx <- order(rtime)
        rtime <- rtime[idx]
        intensity <- intensity[idx]
    }
    new("Chromatogram", rtime = rtime, intensity = intensity,
        mz = range(mz), filterMz = range(filterMz),
        precursorMz = range(precursorMz), productMz = range(productMz),
        fromFile = as.integer(fromFile), aggregationFun = aggregationFun,
        msLevel = as.integer(msLevel))
}

#' @description Plot a single Chromatogram object
#'
#' @author Johannes Rainer
#'
#' @noRd
.plotChromatogram <- function(x, rt, col = "#00000060", lty = 1, type = "l",
                              xlab = "retention time", ylab = "intensity",
                              main = NULL, ...) {
    if (is.null(main)) {
        suppressWarnings(
            mzr <- range(mz(x), na.rm = TRUE, finite = TRUE)
        )
        main <- paste0(format(mzr, digits = 7), collapse = " - ")
    }
    if (!missing(rt))
        x <- filterRt(x, rt = rt)
    plot(x = rtime(x), y = intensity(x), main = main, col = col, lty = lty,
         type = type, xlab = xlab, ylab = ylab, ...)
}

#' @rdname Chromatogram-class
aggregationFun <- function(object) {
    if (!is(object, "Chromatogram"))
        stop("'object' is supposed to be a 'Chromatogram' class")
    object@aggregationFun
}

#' Simple function to bin intensities along retention time for a Chromatogram
#' object.
#'
#' @author Johannes Rainer
#'
#' @noRd
.bin_Chromatogram <- function(x, binSize = 0.5,
                              breaks = seq(floor(min(rtime(x))),
                                           ceiling(max(rtime(x))),
                                           by = binSize),
                              fun = max) {
    bins <- .bin_values(x@intensity, x@rtime, binSize = binSize,
                        breaks = breaks, fun = fun)
    x@intensity <- bins$x
    x@rtime <- bins$mids
    if (validObject(x))
        x
}

.normalize_chromatogram <- function(x, method = "max") {
    ref <- do.call(method, list(x = x@intensity, na.rm = TRUE))
    x@intensity <- x@intensity / ref
    x
}

.filter_intensity_chromatogram <- function(x, intensity = 0, ...) {
    if (is.numeric(intensity)) {
        keep <- x@intensity >= intensity[1]
    } else if (is.function(intensity)) {
        keep <- intensity(x, ...)
    } else stop("'intensity' should be either a numeric value or a function.")
    if (!is.logical(keep) | length(keep) != length(x@intensity))
        stop("The filter function seems to not return the expected result.")
    keep <- which(keep)
    x@intensity <- x@intensity[keep]
    x@rtime <- x@rtime[keep]
    x
}

.align_chromatogram <- function(x, y, method = c("closest", "approx"),
                                na.value = NA_real_, ...) {
    method <- match.arg(method)
    switch(
        method,
        closest = .align_chromatogram_closest(x = x, y = y,
                                              na.value = na.value, ...),
        approx = .align_chromatogram_approx(x = x, y = y,
                                            na.value = na.value, ...))
}

#' @description
#'
#' *Align* chromatogram `x` to chromatogram `y`.
#' The retention times of the first `Chromatogram` (`x`) will be replaced
#' by the retention times of the second (`y`) estimated by the `approx`
#' function.
#'
#' This type of alignment should be used if the chromatograms were measured in
#' the same run (e.g. in SWATH experiments).
#'
#' @param x [Chromatogram()] object that will be aligned against `y`
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value optional parameter allowing to specify the value with which
#'     `NA`s in the aligned chromatogram `x` will be replaced.
#'
#' @return aligned `x` `Chromatogram` with the same number of data points and
#'     same retention times than `y`.
#'
#' @author Michael Witting
#'
#' @noRd
.align_chromatogram_approx <- function(x, y, na.value = NA_real_, ...) {
    x_aligned <- approx(x@rtime, x@intensity, y@rtime)
    if (!is.na(na.value)) {
        x_aligned$y[is.na(x_aligned$y)] <- na.value
    }
    ## correct rtime and int in chromatogram object
    x@rtime <- x_aligned$x
    x@intensity <- x_aligned$y
    x
}

#' @description
#'
#' *Align* chromatogram `x` against chromatogram `y` matching each data point
#' in `x` to the data point in `y` with the smallest difference between their
#' retention times, but only if the difference in retention times is `<=` the
#' minimal average difference between retention times in `x` or `y` (otherwise
#' the data point will be discarded).
#'
#' @param x [Chromatogram()] object that will be aligned against `y`
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value optional parameter allowing to specify the value with which
#'     `NA`s in the aligned chromatogram `x` will be replaced.
#'
#' @param ... optional parameters to be passed along to the `.match_closest`
#'     function.
#'
#' @return aligned `x` `Chromatogram` with the same number of data points and
#'     same retention times than `y`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.align_chromatogram_closest <-
    function(x, y, na.value = NA_real_,
             tolerance = min(mean(diff(x@rtime)), mean(diff(y@rtime))), ...) {
        idx <- closest(x = x@rtime, table = y@rtime, duplicates = "closest",
                       tolerance = tolerance)
        not_na <- !is.na(idx)
        x@rtime <- rtime(y)
        new_int <- rep(na.value, length(y))
        if (any(not_na))
            new_int[idx[not_na]] <- x@intensity[not_na]
        x@intensity <- new_int
        x
    }

.chromatogram_transform_intensity <- function(x, FUN = identity) {
    x@intensity <- FUN(x@intensity)
    x
}
