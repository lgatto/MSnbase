#' @title Representation of chromatographic MS data
#'
#' @aliases compareChromatograms transformIntensity
#'
#' @name Chromatogram-class
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
#' @param ALIGNFUN for `compareChromatograms`: function to align chromatogram
#'     `x` against chromatogram `y`. Defaults to `alignRt`.
#'
#' @param ALIGNFUNARGS `list` of parameters to be passed to `ALIGNFUN`.
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
#' @param FUN for `compareChromatograms`: function to calculate a similarity
#'     score on the intensity values of the compared and aligned chromatograms.
#'     Defaults to `FUN = cor`. For `transformIntensity`: function to transform
#'     chromatograms' intensity values. Defaults to `FUN = identity`.
#'
#' @param FUNARGS for `compareChromatograms`: `list` with additional parameters
#'     for `FUN`. Defaults to `FUNARGS = list(use = "pairwise.complete.obs")`.
#'
#' @param intensity for `Chromatogram`: `numeric` with the intensity values
#'     (length has to be equal to the length of `rtime`). For `filterIntensity`:
#'     `numeric(1)` or `function` to use to filter intensities. See description
#'     for details.
#'
#' @param lty for `plot`: the line type. See help page of `plot` in
#'     the `graphics` package for details.
#'
#' @param main for `plot`: the plot title. If not provided the mz range
#'     will be used as plot title.
#'
#' @param method `character(1)`. For `normalise`: defining whether each
#'     chromatogram should be normalized to its maximum signal
#'     (`method = "max"`) or total signal (`method = "sum"`).
#'     For `alignRt`: aligning approach that should be used (see description).
#'     Defaults to `method = "closest"`.
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
#' @param y for `alignRt`: `Chromatogram` against which `x` should be aligned
#'     against.
#'
#' @param ylab for `plot`: the y-axis label.
#'
#' @param ... for `plot`: additional arguments to be passed to the
#'     base `plot` function. For `filterIntensity`: additional parameters passed
#'     along to the function provided with `intensity`.
#'     For `compareChromatograms`: ignored
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
#' - `filterIntensity`: filter a [Chromatogram()] object removing data
#'   points with intensities below a user provided threshold. If `intensity`
#'   is a `numeric` value, the returned chromatogram will only contain data
#'   points with intensities > `intensity`. In addition it is possible to
#'   provide a function to perform the filtering.
#'   This function is expected to take the input `Chromatogram` (`object`) and
#'   to return a logical vector with the same length then there are data points
#'   in `object` with `TRUE` for data points that should be kept and `FALSE`
#'   for data points that should be removed. See examples below.
#'
#'
#' @section Data processing and manipulation:
#'
#' - `alignRt`: Aligns chromatogram `x` against chromatogram `y`. The resulting
#'   chromatogram has the same length (number of data points) than `y` and the
#'   same retention times thus allowing to perform any pair-wise comparisons
#'   between the chromatograms. If `x` is a [MChromatograms()] object, each
#'   `Chromatogram` in it is aligned against `y`. Additional parameters (`...`)
#'   are passed along to the alignment functions (e.g. [closest()]).
#'
#'   Parameter `method` allows to specify which alignment method
#'   should be used. Currently there are the following options:
#'
#'   - `method = "closest"` (the default): match data points in the first
#'     chromatogram (`x`) to those of the second (`y`) based on the difference
#'     between their retention times: each data point in `x` is assigned to the
#'     data point in `y` with the smallest difference in their retention times
#'     if their difference is smaller than the minimum average difference
#'     between retention times in `x` or `y` (parameter `tolerance` for the
#'     call to the [closest()] function).
#'     By setting `tolerance = 0` only exact retention times are matched against
#'     each other (i.e. only values are kept with exactly the same retention
#'     times between both chromatograms).
#'   - `method = "approx"`: uses the base R `approx` function to approximate
#'     intensities in `x` to the retention times in `y` (using linear
#'     interpolation). This should only be used for chromatograms that were
#'     measured in the same measurement run (e.g. MS1 and corresponding MS2
#'     chromatograms from SWATH experiments).
#'
#' - `bin`: aggregates intensity values from a chromatogram in discrete bins
#'   along the retention time axis and returns a `Chromatogram` object with
#'   the retention time representing the mid-point of the bins and the
#'   intensity the binned signal. Parameters `binSize` and `breaks` allow to
#'   define the binning, `fun` the function which should be used to aggregate
#'   the intensities within a bin.
#'
#' - `compareChromatograms`: calculates a similarity score between 2
#'   chromatograms after aligning them. Parameter `ALIGNFUN` allows to define
#'   a function that can be used to align `x` against `y` (defaults to
#'   `ALIGNFUN = alignRt`). Subsequently, the similarity is calculated on the
#'   aligned intensities with the function provided with parameter `FUN` which
#'   defaults to `cor` (hence by default the Pearson correlation is calculated
#'   between the aligned intensities of the two compared chromatograms).
#'   Additional parameters can be passed to the `ALIGNFUN` and `FUN` with the
#'   parameter `ALIGNFUNARGS` and `FUNARGS`, respectively.
#'
#' - `clean`: removes 0-intensity data points (and `NA` values). See [clean()]
#'   for details.
#'
#' - `normalize`, `normalise`: *normalises* the intensities of a chromatogram by
#'   dividing them either by the maximum intensity (`method = "max"`) or total
#'   intensity (`method = "sum"`) of the chromatogram.
#'
#' - `transformIntensity`: allows to manipulate the intensity values of a
#'   chromatogram using a user provided function. See below for examples.
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
#'
#' ## Data manipulations:
#'
#' ## normalize a chromatogram
#' par(mfrow = c(1, 2))
#' plot(chr)
#' plot(normalize(chr, method = "max"))
#'
#' ## Align chromatograms against each other
#'
#' chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'     intensity = c(3, 5, 14, 30, 24, 6, 2, 1, 1, 0))
#' chr2 <- Chromatogram(rtime = c(2.5, 3.42, 4.5, 5.43, 6.5),
#'     intensity = c(5, 12, 15, 11, 5))
#'
#' plot(chr1, col = "black")
#' points(rtime(chr2), intensity(chr2), col = "blue", type = "l")
#'
#' ## Align chr2 to chr1 without interpolation
#' res <- alignRt(chr2, chr1)
#' rtime(res)
#' intensity(res)
#' points(rtime(res), intensity(res), col = "#00ff0080", type = "l")
#'
#' ## Align chr2 to chr1 with interpolation
#' res <- alignRt(chr2, chr1, method = "approx")
#' points(rtime(res), intensity(res), col = "#ff000080", type = "l")
#' legend("topright", col = c("black", "blue", "#00ff0080","#ff000080"),lty = 1,
#'     legend = c("chr1", "chr2", "chr2 matchRtime", "chr2 approx"))
#'
#'
#' ## Compare Chromatograms. Align chromatograms with `alignRt` and
#' ## method `"approx"`
#' compareChromatograms(chr2, chr1, ALIGNFUNARGS = list(method = "approx"))
#'
#' ## Data filtering
#'
#' chr1 <- Chromatogram(rtime = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
#'     intensity = c(3, 5, 14, 30, 24, 6, 2, 1, 1, 0))
#'
#' ## Remove data points with intensities below 10
#' res <- filterIntensity(chr1, 10)
#' intensity(res)
#'
#' ## Remove data points with an intensity lower than 10% of the maximum
#' ## intensity in the Chromatogram
#' filt_fun <- function(x, prop = 0.1) {
#'     x@intensity >= max(x@intensity, na.rm = TRUE) * prop
#' }
#' res <- filterIntensity(chr1, filt_fun)
#' intensity(res)
#'
#' ## Remove data points with an intensity lower than half of the maximum
#' res <- filterIntensity(chr1, filt_fun, prop = 0.5)
#' intensity(res)
#'
#' ## log2 transform intensity values
#' res <- transformIntensity(chr1, log2)
#' intensity(res)
#' log2(intensity(chr1))
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

#' @rdname Chromatogram-class
setMethod("normalize", "Chromatogram",
          function(object, method = c("max", "sum")) {
              method <- match.arg(method)
              .normalize_chromatogram(object, method)
})

#' @rdname Chromatogram-class
setMethod("filterIntensity", "Chromatogram", function(object,
                                                      intensity = 0, ...) {
    .filter_intensity_chromatogram(object, intensity = intensity, ...)
})

#' @rdname Chromatogram-class
setMethod("alignRt", signature = c(x = "Chromatogram", y = "Chromatogram"),
          function(x, y, method = c("closest", "approx"), ...) {
              .align_chromatogram(x = x, y = y, method = method, ...)
          })

#' @rdname Chromatogram-class
setMethod("compareChromatograms",
          signature = c(x = "Chromatogram", y = "Chromatogram"),
          function(x, y, ALIGNFUN = alignRt, ALIGNFUNARGS = list(),
                   FUN = cor, FUNARGS = list(use = "pairwise.complete.obs"),
                   ...) {
              if (length(x) != length(y) || !all(rtime(x) == rtime(y)))
                  x <- do.call(ALIGNFUN, c(list(x, y), ALIGNFUNARGS))
              do.call(FUN, c(list(x@intensity, y@intensity), FUNARGS))
          })

#' @rdname Chromatogram-class
setMethod("transformIntensity", "Chromatogram", function(object,
                                                         FUN = identity) {
    object <- .chromatogram_transform_intensity(object, FUN = FUN)
    validObject(object)
    object
})
