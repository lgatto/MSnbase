#' @title Representation of chromatographic MS data
#'
#' @aliases Chromatogram-class
#'
#' @description The `Chromatogram` class is designed to store
#'     chromatographic MS data, i.e. pairs of retention time and intensity
#'     values. Instances of the class can be created with the
#'     `Chromatogram` constructor function but in most cases the dedicated
#'     methods for [OnDiskMSnExp-class] and
#'     [MSnExp-class] objects extracting chromatograms should be
#'     used instead (i.e. the [chromatogram()] method).
#'
#' @details The `mz`, `filterMz`, `precursorMz` and
#'     `productMz` are stored as a `numeric(2)` representing a range
#'     even if the chromatogram was generated for only a single ion (i.e. a
#'     single mz value). Using ranges for `mz` values allow this class to
#'     be used also for e.g. total ion chromatograms or base peak chromatograms.
#'
#'     The slots `precursorMz` and `productMz` allow to represent SRM
#'     (single reaction monitoring) and MRM (multiple SRM) chromatograms. As
#'     example, a `Chromatogram` for a SRM transition 273 -> 153 will have
#'     a `@precursorMz = c(273, 273)` and a
#'     `@productMz = c(153, 153)`.
#'
#' @param aggregationFun for `Chromatogram`: `character` string specifying
#'     the function that was used to aggregate intensity values for the same
#'     retention time across the mz range. Supported are `"sum"` (total ion
#'     chromatogram), `"max"` (base peak chromatogram), `"min"` and `"mean"`.
#'
#' @param all for `clean`: `logical(1)` whether all 0 intensities should be
#'     removed. Defaults to `all = FALSE`. See [clean()] for details.
#'
#' @param binSize for `bin`: `numeric(1)` with the size of the bins
#'     (in seconds). Defaults to `binSize = 0.5`.
#'
#' @param breaks for `bin`: `numeric` defining the bins. Usually not
#'     required as the function calculates the bins automatically based on
#'     `binSize`.
#'
#' @param col for `plot`: the color to be used for plotting.
#'
#' @param filter for `mz`: `logical(1)` defining whether the m/z range to filter
#'     the originating object (e.g. `MSnExp` object) should be returned or the
#'     m/z range of the actual data. Defaults to `filter = FALSE`.
#'
#' @param filterMz for `Chromatogram`: `numeric(2)` representing the mz value
#'     range (min, max) that was used to filter the original object on m/z
#'     dimension. If not applicable use `filterMz = c(0, 0)`.
#'
#' @param fromFile for `Chromatogram`: `integer(1)` the index of the file within
#'     the `OnDiskMSnExp` or `MSnExp` from which the chromatogram was extracted.
#'
#' @param fun for `bin`: function to be used to aggregate the intensity
#'     values falling within each bin. Defaults to `fun = max`.
#'
#' @param intensity for `Chromatogram`: `numeric` with the intensity values
#'     (length has to be equal to the length of `rtime`).
#'
#' @param lty for `plot`: the line type. See help page of `plot` in
#'     the `graphics` package for details.
#'
#' @param main for `plot`: the plot title. If not provided the mz range
#'     will be used as plot title.
#'
#' @param msLevel for `Chromatogram`: `integer(1)` with the MS level from
#'     which the chromatogram was extracted.
#'
#' @param mz for `Chromatogram`: `numeric(2)` representing the mz value range
#'     (min, max) on which the chromatogram was created. This is supposed to
#'     contain the *real* range of mz values in contrast to `filterMz`.
#'     If not applicable use `mzrange = c(0, 0)`.
#'
#' @param na.rm for `clean`: if all `NA` intensities should be removed before
#'     cleaning the `Chromatogram`. Defaults to `clean = FALSE`.
#'
#' @param object `Chromatogram` object.
#'
#' @param precursorMz for `Chromatogram`: `numeric(2)` for SRM/MRM transitions.
#'     Represents the mz of the precursor ion. See details for more information.
#'
#' @param productMz for `Chromatogram`: `numeric(2)` for SRM/MRM transitions.
#'     Represents the mz of the product. See details for more information.
#'
#' @param rt for `filterRt`: `numeric(2)` defining the lower and upper retention
#'     time to which the `Chromatogram` should be subsetted.
#'
#' @param rtime for `Chromatogram`: `numeric` with the retention times (length
#'     has to be equal to the length of `intensity`).
#'
#' @param type for `plot`: the type of plot. See help page of `plot` in
#'     the `graphics` package for details.
#'
#' @param x `Chromatogram` object.
#'
#' @param xlab for `plot`: the x-axis label.
#'
#' @param ylab for `plot`: the y-axis label.
#'
#' @param ... for `plot`: additional arguments to be passed to the
#'     base `plot` function.
#'
#'
#' @section Object creation:
#'
#' `Chromatogram` objects can be extracted from an `MSnExp` or `OnDiskMSnExp`
#' object with the `chromatogram()` function.
#'
#' Alternatively, the constructor function `Chromatogram` can be used, which
#' takes arguments `rtime`, `intensity`, `mz`, `filterMz`, `precursorMz`,
#' `productMz`, `fromFile`, `aggregationFun` and `msLevel`.
#'
#'
#' @section Data access and coercion:
#'
#' - `aggregationFun`: gets the aggregation function used to create the
#'   `Chromatogram`.
#'
#' - `as.data.frame`: returns a `data.frame` with columns `"rtime"` and
#'   `"intensity"`.
#'
#' - `fromFile`: returns an `integer(1)` with the index of the originating file.
#'
#' - `intensity`: returns the intensities from the `Chromatogram`.
#'
#' - `isEmpty`: returns `TRUE` if the chromatogram is empty or has only `NA`
#'    intensities.
#'
#' - `length`: returns the length (i.e. number of data points) of the
#'   `Chromatogram`.
#'
#' - `msLevel`: returns an `integer(1)` with the MS level of the chromatogram.
#'
#' - `mz`: get the m/z (range) from the `Chromatogram`. The function returns
#'   a `numeric(2)` with the lower and upper boundaries. Parameter `filter`
#'   allows to specify whether the m/z range used to filter the originating
#'   object should be returned or the m/z range of the actual data.
#'
#' - `precursorMz`: get the m/z of the precursor ion. The function returns a
#'   `numeric(2)` with the lower and upper boundary.
#'
#' - `productMz`: get the m/z of the producto chromatogram/ion. The function
#'   returns a `numeric(2)` with the lower and upper m/z value.
#'
#' - `rtime`: returns the retention times from the `Chromatogram`.
#'
#'
#' @section Data subsetting and filtering:
#'
#' - `filterRt`: filter/subset the `Chromatogram` to the specified retention
#'   time range (defined with parameter `rt`).
#'
#'
#' @section Data processing and manipulation:
#'
#' - `bin`: aggregates intensity values from a chromatogram in discrete bins
#'   along the retention time axis and returns a `Chromatogram` object with
#'   the retention time representing the mid-point of the bins and the
#'   intensity the binned signal. Parameters `binSize` and `breaks` allow to
#'   define the binning, `fun` the function which should be used to aggregate
#'   the intensities within a bin.
#'
#' - `clean`: removes 0-intensity data points (and `NA` values). See [clean()]
#'   for details.
#'
#'
#' @section Data visualization:
#'
#' - `plot`: plots a `Chromatogram` object.
#'
#'
#' @rdname Chromatogram-class
#'
#' @export
#'
#' @seealso [MChromatograms] for combining `Chromatogram` in
#'     a two-dimensional matrix (rows being mz-rt ranges, columns samples).
#'     `chromatogram()] for the method to extract chromatogram data
#'     from an `MSnExp` or `OnDiskMSnExp` object.
#'     [clean()] for the method to *clean* a `Chromatogram` object.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Create a simple Chromatogram object.
#' ints <- abs(rnorm(100, sd = 100))
#' rts <- seq_len(length(ints))
#' chr <- Chromatogram(rtime = rts, intensity = ints)
#' chr
#'
#' ## Extract intensities
#' intensity(chr)
#'
#' ## Extract retention times
#' rtime(chr)
#'
#' ## Extract the mz range - is NA for the present example
#' mz(chr)
#'
#' ## plot the Chromatogram
#' plot(chr)
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
NULL

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

#' @rdname Chromatogram-class
setMethod("rtime", "Chromatogram", function(object) {
    object@rtime
})

#' @rdname Chromatogram-class
setMethod("intensity", "Chromatogram", function(object) {
    object@intensity
})

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

#' @rdname Chromatogram-class
setMethod("precursorMz", "Chromatogram", function(object) {
    object@precursorMz
})

#' @rdname Chromatogram-class
setMethod("fromFile", "Chromatogram", function(object) {
    object@fromFile
})

#' @rdname Chromatogram-class
setMethod("length", "Chromatogram", function(x) {
    length(x@rtime)
})

#' @rdname Chromatogram-class
setMethod("as.data.frame", "Chromatogram", function(x) {
    data.frame(rtime = x@rtime, intensity = x@intensity)
})

#' @rdname Chromatogram-class
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

#' @rdname Chromatogram-class
setMethod("msLevel", "Chromatogram", function(object) {
    object@msLevel
})

#' @rdname Chromatogram-class
setMethod("isEmpty", "Chromatogram", function(x) {
    (length(x) == 0 | all(is.na(intensity(x))))
})

#' @rdname Chromatogram-class
setMethod("productMz", "Chromatogram", function(object) {
    object@productMz
})

#' @rdname Chromatogram-class
setMethod("bin", "Chromatogram", .bin_Chromatogram)
