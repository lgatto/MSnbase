#' @include functions-MSpectra.R

setMethod("show", "MSpectra", function(object) {
    .show_MSpectra(object, margin = "  ", print.classinfo = TRUE)
})

#' @rdname MSpectra
#'
#' @section Accessing spectrum attributes:
#'
#' These methods allow to access the attributes and values of the individual
#' `Spectrum` ([Spectrum1-class] or [Spectrum2-class]) objects within the list.
#'
#' - `mz` return the m/z values of each spectrum as a `list` of `numeric`
#'   vectors.
#'
#' - `intensity` return the intensity values of each spectrum as a `list` of
#'   `numeric` vectors.
#'
#' - `rtime` return the retention time of each spectrum as a `numeric` vector
#'   with length equal to the length of `object`.
#'
#' - `precursorMz`, `precursorCharge`, `precursorIntensity`, `precScanNum`
#'   return precursor m/z values, charge, intensity and scan number for each
#'   spectrum as a `numeric` (or `integer`) vector with length equal to the
#'   length of `object`. Note that for [Spectrum1-class] objects `NA` will be
#'   returned.
#'
#' - `acquisitionNum` and `scanIndex` return the acquisition number of each
#'   spectrum and its scan index as an `integer` vector with the same length
#'   than `object`.
#'
#' - `ionCount` and `tic` return the ion count and total ion current of each
#'   spectrum.
#'
#' - `peaksCount` returns the number of peaks for each spectrum as a `integer`
#'   vector.
#'
#' - `msLevel` returns the MS level of each spectrum.
#'
#' - `collisionEnergy` returns the collision energy for each spectrum or `NA`
#'   for [Spectrum1-class] objects.
#'
#' - `polarity` returns the spectra's polarity.
#'
#' - `fromFile` returns the index from the (e.g. mzML) file the spectra where
#'   from. This applies only for spectra read using the [readMSData()] function.
#'
#' - `smoothed` whether spectra have been smoothed (i.e. processed with the
#'   [smooth()] method. Returns a `logical` of length equal to the
#'   number of spectra.
#'
#' - `isEmpty` returns `TRUE` for spectra without peak data.
#'
#' - `centroided`, `isCentroided` returns for each spectrum whether it contains
#'   *centroided* data. While `centroided` returns the internal attribute of
#'   each spectrum, `isCentroided` tries to guess whether spectra are
#'   centroided from the actual peak data.
#'
#' @section Data manipulation methods:
#'
#' - `clean` *cleans* each spectrum. See [clean()] for more details.
#'
#' - `pickPeaks` performs peak picking to generate centroided spectra. See
#'   [pickPeaks()] for more details.
#'
#' - `removePeaks` removes peaks lower than a threshold `t`. See
#'   [removePeaks()] for more details.
#'
#' - `smooth` *smooths* spectra. See [smooth()] for more details.
#'
#' @section Filtering and subsetting:
#'
#' - `[` can be used to subset the `MSpectra` object.
#'
#' - `filterMsLevel` filters `MSpectra` to retain only spectra from certain MS
#'   level(s).
#'
#' - `filterMz` filters the spectra by the specified `mz` range. See
#'   [filterMz()] for details.
#'
#'
#' @md
#'
#' @param all For `clean`: if `FALSE` original 0-intensity values are retained
#'   around peaks.
#'
#' @param msLevel. For `clean`, `removePeaks`, `filterMz`, `pickPeaks`:
#'   optionally specify the MS level(s) of the spectra on which the operation
#'   should be performed.
#'   For `filterMsLevels`: MS level(s) to which the `MSpectra` should be reduced.
#'
#' @param method For `pickPeaks` and `smooth`: see [pickPeaks()] and [smooth()]
#'   for details.
#'
#' @param mz For `filterMz`: `numeric(2)` defining the lower and upper m/z
#'   for the filter. See [filterMz()] for details.
#'
#' @param t For `removePeaks`: `numeric(1)` specifying the threshold below
#'   which intensities are set to 0.
#'
#' @param halfWindowSize For `pickPeaks` and `smooth`: see [pickPeaks()]
#'   and [smooth()] for details.
#'
#' @param SNR For `pickPeaks`: see [pickPeaks()] for details.
#'
#' @param refineMz For `pickPeaks`: see [pickPeaks()] for details.
#'
#' @examples
#'
#' ## Extract the mz values for the individual spectra
#' mz(spl)
setMethod("mz", "MSpectra", function(object) {
    lapply(object, function(z) z@mz)
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the intensity values for the individual spectra
#' intensity(spl)
setMethod("intensity", "MSpectra", function(object) {
    lapply(object, function(z) z@intensity)
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the retention time values for the individual spectra
#' rtime(spl)
setMethod("rtime", "MSpectra", function(object) {
    vapply(object, function(z) if(length(z@rt)) z@rt else NA_real_, numeric(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the precursor m/z of each spectrum.
#' precursorMz(spl)
setMethod("precursorMz", "MSpectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precursorMz)) z@precursorMz
        else NA_real_
    }, numeric(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the precursor charge of each spectrum.
#' precursorCharge(spl)
setMethod("precursorCharge", "MSpectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precursorCharge)) z@precursorCharge
        else NA_integer_
    }, integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the precursor scan number for each spectrum.
#' precScanNum(spl)
setMethod("precScanNum", "MSpectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precScanNum)) z@precScanNum
        else NA_integer_
    }, integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the precursor intensity of each spectrum.
#' precursorIntensity(spl)
setMethod("precursorIntensity", "MSpectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@precursorIntensity))
            z@precursorIntensity
        else NA_real_
    }, numeric(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the acquisition number of each spectrum.
#' acquisitionNum(spl)
setMethod("acquisitionNum", "MSpectra", function(object) {
    vapply(object, function(z)
        if (length(z@acquisitionNum)) z@acquisitionNum else NA_integer_,
        integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the scan index of each spectrum.
#' scanIndex(spl)
setMethod("scanIndex", "MSpectra", function(object) {
    vapply(object, function(z)
        if (length(z@scanIndex)) z@scanIndex else NA_integer_,
        integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Get the number of peaks per spectrum.
#' peaksCount(spl)
setMethod("peaksCount", "MSpectra", function(object) {
    vapply(object, peaksCount, integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Get the MS level of each spectrum.
#' msLevel(spl)
setMethod("msLevel", "MSpectra", function(object) {
    vapply(object, msLevel, integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Get the total ion current for each spectrum.
#' tic(spl)
setMethod("tic", "MSpectra", function(object) {
    vapply(object, tic, numeric(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Get the total ion current for each spectrum.
#' ionCount(spl)
setMethod("ionCount", "MSpectra", function(object) {
    vapply(object, ionCount, numeric(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the collision energy for each spectrum.
#' collisionEnergy(spl)
setMethod("collisionEnergy", "MSpectra", function(object) {
    vapply(object, function(z) {
        if (is(z, "Spectrum2") && length(z@collisionEnergy)) z@collisionEnergy
        else NA_real_
    }, numeric(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Extract the file index for each spectrum.
#' fromFile(spl)
setMethod("fromFile", "MSpectra", function(object) {
    vapply(object, function(z)
        if (length(z@fromFile)) z@fromFile else NA_integer_,
        integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Get the polarity for each spectrum.
#' polarity(spl)
setMethod("polarity", "MSpectra", function(object) {
    vapply(object, function(z)
        if (length(z@polarity)) z@polarity else NA_integer_,
        integer(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Whether spectra are smoothed (i.e. processed with the `smooth`
#' ## function).
#' smoothed(spl)
setMethod("smoothed", "MSpectra", function(object) {
    vapply(object, smoothed, logical(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Are spectra empty (i.e. contain no peak data)?
#' isEmpty(spl)
setMethod("isEmpty", "MSpectra", function(x) {
    vapply(x, isEmpty, logical(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Do the spectra contain centroided data?
#' centroided(spl)
setMethod("centroided", "MSpectra", function(object) {
    vapply(object, centroided, logical(1))
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Do the spectra contain centroided data? Whether spectra are centroided
#' ## is estimated from the peak data.
#' isCentroided(spl)
setMethod("isCentroided", "MSpectra", function(object) {
    vapply(object, isCentroided, logical(1))
})

#' @rdname MSpectra
#'
#' @description
#'
#' `writeMgfData` exports a `MSpectra` object to a file in MGF format. All
#' metadata columns present in `mcols` are exported as additional fields with
#' the capitalized column names used as field names (see examples below).
#'
#' @param con For `writeMgfData`: `character(1)` defining the file name of
#'     the MGF file.
#'
#' @param COM For `writeMgfData`: optional `character(1)` providing a comment
#'     to be added to the file.
#'
#' @param TITLE For `writeMgfData`: optional `character(1)` defining the title
#'     for the MGF file.
#'
#' @md
#'
#' @examples
#'
#' ## Export the spectrum list to a MGF file. Values in metadata columns are
#' ## exported as additional field for each spectrum.
#' tmpf <- tempfile()
#' writeMgfData(spl, tmpf)
#'
#' ## Evaluate the written output. The ID of each spectrum (defined in the
#' ## "id" metadata column) is exported as field "ID".
#' readLines(tmpf)
#'
#' ## Set mcols to NULL to avoid export of additional data fields.
#' mcols(spl) <- NULL
#' file.remove(tmpf)
#'
#' writeMgfData(spl, tmpf)
#' readLines(tmpf)
setMethod("writeMgfData", "MSpectra", function(object, con = "spectra.mgf",
                                              COM = NULL, TITLE = NULL) {
    if (file.exists(con))
        stop("file ", con, " does already exist.")
    writeMgfDataFile(as.list(object), con = con, COM = COM,
                     TITLE = TITLE, addFields = mcols(object))
})

#' @rdname MSpectra
setMethod("clean", "MSpectra", function(object, all = FALSE,
                                       msLevel. = msLevel., ...) {
    object <- endoapply(object, clean, all = all, msLevel. = msLevel., ...)
    if (validObject(object))
        object
})

#' @rdname MSpectra
setMethod("removePeaks", "MSpectra", function(object, t, msLevel., ...) {
    object <- endoapply(object, removePeaks, t = t, msLevel. = msLevel., ...)
    if (validObject(object))
        object
})

#' @rdname MSpectra
setMethod("filterMz", "MSpectra", function(object, mz, msLevel., ...) {
    object <- endoapply(object, filterMz, mz = mz, msLevel. = msLevel., ...)
    if (validObject(object))
        object
})

#' @rdname MSpectra
setMethod("pickPeaks", "MSpectra", function(object, halfWindowSize = 3L,
                                           method = c("MAD", "SuperSmoother"),
                                           SNR = 0L,
                                           refineMz = c("none", "kNeighbors",
                                                        "kNeighbours",
                                                        "descendPeak"),
                                           msLevel. = unique(msLevel(object)),
                                           ...) {
    object <- endoapply(object, pickPeaks, halfWindowSize = halfWindowSize,
                        method = match.arg(method), SNR = SNR,
                        refineMz = refineMz, msLevel. = msLevel., ...)
    if (validObject(object))
        object
})

#' @rdname MSpectra
setMethod("smooth", "MSpectra", function(x, method = c("SavitzkyGolay",
                                                      "MovingAverage"),
                                        halfWindowSize = 2L, ...) {
    x <- endoapply(x, smooth, method = match.arg(method),
                   halfWindowSize = halfWindowSize, ...)
    if (validObject(x))
        x
})

#' @title Combine Spectra
#'
#' @description
#'
#' `combineSpectra` combines spectra in a [MSnExp-class], [OnDiskMSnExp-class]
#' or [MSpectra-class] object applying the summarization function `fun` to sets
#' of spectra defined by a factor (`fcol` parameter). The resulting combined
#' spectrum for each set contains metadata information (present in `mcols` and
#' all spectrum information other than `mz` and `intensity`) from the **first**
#' spectrum in each set.
#'
#' Combining of spectra for [MSnExp-class] or [OnDiskMSnExp-class] objects is
#' performed by default for each file **separately**, combining of spectra
#' across files is thus not possible. See examples for details.
#'
#' @aliases combineSpectra
#'
#' @rdname combineSpectra
#'
#' @param object A [MSnExp-class] or [MSpectra-class]
#'
#' @param fcol For `MSpectra` objects: `mcols` column name to be used to define
#'     the sets of spectra to be combined. If missing, all spectra are
#'     considered to be one set. For `MSnExp`/`OnDiskMSnExp` objects: column
#'     in `fData(object)` defining which spectra to combine. See examples below
#'     for more details.
#'
#' @param fun *Deprecated* use `method` instead.
#'
#' @param method `function` to be used to combine the spectra by `fcol`. Has to
#'     be a function that takes a list of spectra as input and returns a single
#'     [Spectrum-class]. See [meanMzInts()] for details.
#'
#' @param ... additional arguments for `fun`.
#'
#' @param BPPARAM For `MSnExp`/`OnDiskMSnExp` objects: parallel processing setup
#'     to perform per-file parallel spectra combining. See [bpparam()] for more
#'     details.
#'
#' @return A `MSpectra` or `MSnExp` object with combined spectra. Metadata
#'     (`mcols`) and all spectrum attributes other than `mz` and `intensity`
#'     are taken from the first `Spectrum` in each set.
#'
#' @md
#'
#' @author Johannes Rainer, Laurent Gatto
#'
#' @seealso [meanMzInts()] for a function to combine spectra.
#'
#' @examples
#'
#' set.seed(123)
#' mzs <- seq(1, 20, 0.1)
#' ints1 <- abs(rnorm(length(mzs), 10))
#' ints1[11:20] <- c(15, 30, 90, 200, 500, 300, 100, 70, 40, 20) # add peak
#' ints2 <- abs(rnorm(length(mzs), 10))
#' ints2[11:20] <- c(15, 30, 60, 120, 300, 200, 90, 60, 30, 23)
#' ints3 <- abs(rnorm(length(mzs), 10))
#' ints3[11:20] <- c(13, 20, 50, 100, 200, 100, 80, 40, 30, 20)
#'
#' ## Create the spectra.
#' sp1 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
#'     intensity = ints1, rt = 1)
#' sp2 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
#'     intensity = ints2, rt = 2)
#' sp3 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.009),
#'     intensity = ints3, rt = 3)
#'
#' spctra <- MSpectra(sp1, sp2, sp3,
#'     elementMetadata = DataFrame(idx = 1:3, group = c("b", "a", "a")))
#'
#' ## Combine the spectra reporting the maximym signal
#' res <- combineSpectra(spctra, mzd = 0.05, intensityFun = max)
#' res
#'
#' ## All values other than m/z and intensity are kept from the first spectrum
#' rtime(res)
#'
#' ## Plot the individual and the merged spectrum
#' par(mfrow = c(2, 1), mar = c(4.3, 4, 1, 1))
#' plot(mz(sp1), intensity(sp1), xlim = range(mzs[5:25]), type = "h", col = "red")
#' points(mz(sp2), intensity(sp2), type = "h", col = "green")
#' points(mz(sp3), intensity(sp3), type = "h", col = "blue")
#' plot(mz(res[[1]]), intensity(res[[1]]), type = "h",
#'     col = "black", xlim = range(mzs[5:25]))
#'
#' ## Combine spectra in two sets.
#' res <- combineSpectra(spctra, fcol = "group", mzd = 0.05)
#' res
#'
#' rtime(res)
#'
#' ## Plot the individual and the merged spectra
#' par(mfrow = c(3, 1), mar = c(4.3, 4, 1, 1))
#' plot(mz(sp1), intensity(sp1), xlim = range(mzs[5:25]), type = "h", col = "red")
#' points(mz(sp2), intensity(sp2), type = "h", col = "green")
#' points(mz(sp3), intensity(sp3), type = "h", col = "blue")
#' plot(mz(res[[1]]), intensity(res[[1]]), xlim = range(mzs[5:25]), type = "h",
#'     col = "black")
#' plot(mz(res[[2]]), intensity(res[[2]]), xlim = range(mzs[5:25]), type = "h",
#'     col = "black")
#'
#' ## Combining spectra of an MSnExp/OnDiskMSnExp objects
#' ## Reading data from 2 mzML files
#' sciex <- readMSData(dir(system.file("sciex", package = "msdata"),
#'     full.names = TRUE), mode = "onDisk")
#'
#' ## Filter the file to a retention time range from 2 to 20 seconds (to reduce
#' ## execution time of the example)
#' sciex <- filterRt(sciex, rt = c(2, 20))
#' table(fromFile(sciex))
#'
#' ## We have thus 64 spectra per file.
#'
#' ## In the example below we combine spectra measured in one second to a
#' ## single spectrum. We thus first define the grouping variable and add that
#' ## to the `fData` of the object. For combining, we use the
#' ## `consensusSpectrum` function that combines the spectra keeping only peaks
#' ## that were found in 50% of the spectra; by defining `mzd = 0.01` all peaks
#' ## within an m/z of 0.01 are evaluated for combining.
#' seconds <- round(rtime(sciex))
#' head(seconds)
#' fData(sciex)$second <- seconds
#'
#' res <- combineSpectra(sciex, fcol = "second", mzd = 0.01, minProp = 0.1,
#'     method = consensusSpectrum)
#' table(fromFile(res))
#'
#' ## The data was reduced to 19 spectra for each file.
setMethod("combineSpectra", "MSpectra", function(object, fcol,
                                                 method = meanMzInts, fun, ...) {
    if (!missing(fun)) {
        .Deprecated(
            msg = "Parameter 'fun' is deprecated. Please use 'method' instead")
        if (missing(method))
            method <- fun
    }
    if (missing(fcol)) {
        .by <- factor(rep(1, length(object)))
    } else {
        if (!any(fcol %in% colnames(mcols(object))))
            stop("'fcol' does not match any column names of 'mcols(object)'")
        .by <- factor(mcols(object)[, fcol],
                      levels = unique(mcols(object)[, fcol]))
    }
    res <- lapply(split(object, .by), FUN = method, ...)
    elm <-  mcols(object, use.names = TRUE)[
        levelIndex(.by, which = "first"), , drop = FALSE]
    names(res) <- rownames(elm)
    MSpectra(res, elementMetadata = elm)
})

setAs("MSpectra", "list", function(from) {
    from@listData
})

setAs("MSpectra", "MSnExp", function(from) {
    if (length(unique(msLevel(from))) > 1)
        stop("'from' contains Spectra from more than one MS level. Use ",
             "'filterMsLevel' to restrict to Spectra from a single MS level ",
             "before coercing.")
    if (!length(names(from)))
        names(from) <- as.character(seq_len(length(from)))
    if (is.null(mcols(from)))
        fd <- AnnotatedDataFrame(data.frame(spectrum = seq_len(length(from)),
                                            row.names = names(from)))
    else
        fd <- AnnotatedDataFrame(as.data.frame(mcols(from)))
    fromf <- as.character(unique(fromFile(from)))
    process <- new("MSnProcess", files = fromf,
                   processing = paste("Data converted from Spectra:", date()))
    pd <- data.frame(sampleNames = fromf)
    assaydata <- list2env(as(from, "list"))
    lockEnvironment(assaydata, bindings = TRUE)
    new("MSnExp", assayData = assaydata,
        phenoData = new("AnnotatedDataFrame", pd), featureData = fd,
        processingData = process)
})

#' @rdname MSpectra
#'
#' @examples
#'
#' ## Filter the object by MS level
#' filterMsLevel(spl, msLevel. = 1)
setMethod("filterMsLevel", "MSpectra", function(object, msLevel.) {
    if (missing(msLevel.)) return(object)
    msLevel. <- as.numeric(msLevel.)
    object[msLevel(object) %in% msLevel.]
})

## Still to implement:
## quantify, method = c("trapezoidation", "max", "sum", reporters, strict = FALSE)
## normalize, method = c("max", "sum", "precursor", precursorIntensity)
## bin, binSize = 1L, breaks....
