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

#' @description \code{Chromatogram}: create an instance of the
#'     \code{Chromatogram} class.
#'
#' @param rtime \code{numeric} with the retention times (length has to be equal
#'     to the length of \code{intensity}).
#'
#' @param intensity \code{numeric} with the intensity values (length has to be
#'     equal to the length of \code{rtime}).
#'
#' @param mz \code{numeric(2)} representing the mz value range (min, max)
#'     on which the chromatogram was created. This is supposed to contain the
#'     \emph{real} range of mz values in contrast to the \code{filterMz} below.
#'     If not applicable use \code{mzrange = c(0, 0)}.
#'
#' @param filterMz \code{numeric(2)} representing the mz value range (min,
#'     max) that was used to filter the original object on mz dimension. If not
#'     applicable use \code{filterMz = c(0, 0)}.
#'
#' @param precursorMz \code{numeric(2)} for SRM/MRM transitions.
#'     Represents the mz of the precursor ion. See details for more information.
#' 
#' @param productMz \code{numeric(2)} for SRM/MRM transitions.
#'     Represents the mz of the product. See details for more information.
#' 
#' @param fromFile \code{integer(1)} the index of the file within the
#'     \code{\linkS4class{OnDiskMSnExp}} or \code{\linkS4class{MSnExp}}
#'     from which the chromatogram was extracted.
#'
#' @param aggregationFun \code{character} string specifying the function that
#'     was used to aggregate intensity values for the same retention time across
#'     the mz range. Supported are \code{"sum"} (total ion chromatogram),
#'     \code{"max"} (base peak chromatogram), \code{"min"} and \code{"mean"}.
#'
#' @param msLevel \code{integer} with the MS level from which the chromatogram
#'     was extracted.
#' 
#' @slot .__classVersion__,rtime,intensity,mz,filterMz,precursorMz,productMz,fromFile,aggregationFun,msLevel See corresponding parameter above.
#' 
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

#' @aliases aggregationFun
#'
#' @description \code{aggregationFun,aggregationFun<-} get or set the
#'     aggregation function.
#' 
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
.bin_Chromatogram <- function(object, binSize = 0.5,
                              breaks = seq(floor(min(rtime(object))),
                                           ceiling(max(rtime(object))),
                                           by = binSize),
                              fun = max) {
    bins <- .bin_values(object@intensity, object@rtime, binSize = binSize,
                        breaks = breaks, fun = fun)
    object@intensity <- bins$x
    object@rtime <- bins$mids
    if (validObject(object))
        object
}

