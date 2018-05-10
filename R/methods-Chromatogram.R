setMethod("initialize", "Chromatogram", function(.Object, ...) {
    classVersion(.Object)["Chromatogram"] <- "0.0.1"
    callNextMethod(.Object, ...)
})


#' @rdname Chromatogram-class
setMethod("show", "Chromatogram", function(object) {
    cat("Object of class: ", class(object), "\n", sep = "")
    if (length(object@aggregationFun))
        cat("Intensity values aggregated using:", object@aggregationFun, "\n")
    ## if (length(object@aggregationFun))
    ##     cat(names(.SUPPORTED_AGG_FUN_CHROM)[.SUPPORTED_AGG_FUN_CHROM ==
    ##                                         object@aggregationFun], "\n")
    cat("length of object: ", length(object@rtime), "\n", sep = "")
    cat("from file: ", object@fromFile, "\n", sep = "")
    cat("mz range: [", object@mz[1], ", ", object@mz[2], "]\n", sep = "")
    if (length(object@rtime) > 0) {
        rtr <- range(object@rtime)
        cat("rt range: [", rtr[1], ", ", rtr[2], "]\n", sep = "")
    }
    cat("MS level: ", paste(object@msLevel, collapse = ", "), "\n", sep = "")
})

#' @description \code{rtime} returns the retention times for the rentention time
#'     - intensity pairs stored in the chromatogram.
#'
#' @param object A \code{Chromatogram} object.
#' 
#' @rdname Chromatogram-class
setMethod("rtime", "Chromatogram", function(object) {
    object@rtime
})

#' @description \code{intensity} returns the intensity for the rentention time
#'     - intensity pairs stored in the chromatogram.
#' 
#' @rdname Chromatogram-class
setMethod("intensity", "Chromatogram", function(object) {
    object@intensity
})

#' @description \code{mz} get the mz (range) of the chromatogram. The
#'     function returns a \code{numeric(2)} with the lower and upper mz value.
#'
#' @param filter For \code{mz}: whether the mz range used to filter the
#'     original object should be returned (\code{filter = TRUE}), or the mz
#'     range calculated on the real data (\code{filter = FALSE}).
#' 
#' @rdname Chromatogram-class
setMethod("mz", "Chromatogram", function(object, filter = FALSE) {
    if (filter)
        return(object@filterMz)
    object@mz
})
## #' @rdname Chromatogram-class
## setReplaceMethod("mz", "CentWaveParam", function(object, value) {
##     object@mzrange <- value
##     if (validObject(object))
##         return(object)
## })

#' @description \code{precursorMz} get the mz of the precursor ion. The
#'     function returns a \code{numeric(2)} with the lower and upper mz value.
#' 
#' @rdname Chromatogram-class
setMethod("precursorMz", "Chromatogram", function(object) {
    object@precursorMz
})

#' @description \code{fromFile} returns the value from the \code{fromFile} slot.
#' 
#' @rdname Chromatogram-class
setMethod("fromFile", "Chromatogram", function(object) {
    object@fromFile
})

#' @description \code{length} returns the length (number of retention time -
#'     intensity pairs) of the chromatogram.
#' 
#' @param x For \code{as.data.frame} and \code{length}: a \code{Chromatogram}
#'     object.
#'
#' @rdname Chromatogram-class
setMethod("length", "Chromatogram", function(x) {
    length(x@rtime)
})

#' @description \code{as.data.frame} returns the \code{rtime} and
#'     \code{intensity} values from the object as \code{data.frame}.
#' 
#' @rdname Chromatogram-class
setMethod("as.data.frame", "Chromatogram", function(x) {
    data.frame(rtime = x@rtime, intensity = x@intensity)
})

#' @description \code{filterRt}: filters the chromatogram based on the provided
#'     retention time range.
#'
#' @param rt For \code{filterRt}: \code{numeric(2)} defining the lower and
#'     upper retention time for the filtering.
#'
#' @rdname Chromatogram-class
#'
#' @examples
#'
#' ## Create a simple Chromatogram object based on random values.
#' chr <- Chromatogram(intensity = abs(rnorm(1000, mean = 2000, sd = 200)),
#'         rtime = sort(abs(rnorm(1000, mean = 10, sd = 5))))
#' chr
#'
#' ## Get the intensities
#' head(intensity(chr))
#'
#' ## Get the retention time
#' head(rtime(chr))
#'
#' ## What is the retention time range of the object?
#' range(rtime(chr))
#'
#' ## Filter the chromatogram to keep only values between 4 and 10 seconds
#' chr2 <- filterRt(chr, rt = c(4, 10))
#'
#' range(rtime(chr2))
setMethod("filterRt", "Chromatogram", function(object, rt) {
    if (missing(rt))
        return(object)
    rt <- range(rt)
    ## Use which to be robust against NAs
    keep_em <- which(rtime(object) >= rt[1] & rtime(object) <= rt[2])
    if (length(keep_em)) {
        object@rtime <- rtime(object)[keep_em]
        object@intensity <- intensity(object)[keep_em]
    } else {
        object@rtime <- numeric()
        object@intensity <- numeric()
    }
    if (validObject(object))
        object
})

#' @description \code{clean}: Removes unused 0-intensity data points. See
#'     \code{\link{clean}} documentation for more details and examples.
#'
#' @param all For \code{clean}: \code{logical(1)} whether all 0 intensities
#'     should be removed (default is \code{FALSE}). See \code{\link{clean}} for
#'     more details and examples.
#'
#' @param na.rm For \code{clean}: \code{logical(1)} whether all \code{NA}
#'     intensities should be removed before cleaning the \code{Chromatogram}.
#'     Defaults to \code{FALSE}. See \code{\link{clean}} for more details and
#'     examples.
#' 
#' @rdname Chromatogram-class
setMethod("clean", signature = signature("Chromatogram"),
          function(object, all = FALSE, na.rm = FALSE) {
              keep <- utils.clean(object@intensity, all, na.rm = na.rm)
              object@intensity <- object@intensity[keep]
              object@rtime <- object@rtime[keep]
              if (validObject(object))
                  object
          })

#' @rdname Chromatogram-class
#' 
#' @description \code{plot}: plots a \code{Chromatogram} object.
#'
#' @param col For \code{plot}: the color to be used for plotting.
#'
#' @param lty For \code{plot}: the line type. See \code{\link[graphics]{plot}}
#'     for more details.
#'
#' @param type For \code{plot}: the type of plot. See
#'     \code{\link[graphics]{plot}} for more details.
#'
#' @param xlab For \code{plot}: the x-axis label.
#'
#' @param ylab For \code{plot}: the y-axis label.
#'
#' @param main For \code{plot}: the plot title. If not provided the mz range
#'     will be used as plot title.
#'
#' @param ... For \code{plot}: additional arguments to be passed to the
#'     \code{\link[graphics]{plot}} function.
setMethod("plot", signature = signature("Chromatogram"),
          function(x, col = "#00000060", lty = 1, type = "l",
                   xlab = "retention time", ylab = "intensity",
                   main = NULL, ...) {
              if (isEmpty(x)) {
                  ## Show a warning and plot an empty plot (issue #249)
                  warning("Chromatogram is empty")
                  plot(3, 3, xlab = xlab, ylab = ylab, main = main, pch = NA)
                  text(3, 3, labels = "Empty Chromatogram", col = "red")
              } else {
                  .plotChromatogram(x = x, col = col, lty = lty, type = type,
                                    xlab = xlab, ylab = ylab, main = main, ...)
              }
          })

#' @description \code{msLevel} returns the MS level of the chromatogram.
#' 
#' @rdname Chromatogram-class
setMethod("msLevel", "Chromatogram", function(object) {
    object@msLevel
})

#' @description \code{isEmpty} returns \code{TRUE} for empty chromatogram or
#'     chromatograms with all intensities being \code{NA}.
#' 
#' @rdname Chromatogram-class
setMethod("isEmpty", "Chromatogram", function(x) {
    (length(x) == 0 | all(is.na(intensity(x))))
})

#' @aliases productMz
#' 
#' @description \code{productMz} get the mz of the product chromatogram/ion. The
#'     function returns a \code{numeric(2)} with the lower and upper mz value.
#' 
#' @rdname Chromatogram-class
setMethod("productMz", "Chromatogram", function(object) {
    object@productMz
})

#' @description
#'
#' \code{bin} aggregates intensity values from a chromatogram in discrete bins
#' along the retention time axis and returns a \code{Chromatogram} object with
#' the retention time representing the mid-point of the bins and the intensity
#' the binned signal.
#'
#' @param binSize for \code{bin}: \code{numeric(1)} with the size of the bins
#'     (in seconds).
#'
#' @param breaks for \code{bin}: \code{numeric} defining the bins. Usually not
#'     required as the function calculates the bins automatically based on
#'     \code{binSize}.
#'
#' @param fun for \code{bin}: function to be used to aggregate the intensity
#'     values falling within each bin.
#'
#' @rdname Chromatogram-class
setMethod("bin", "Chromatogram", .bin_Chromatogram)
