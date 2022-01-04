#' @title Container for multiple Chromatogram objects
#'
#' @name MChromatograms-class
#'
#' @aliases coerce,matrix,MChromatograms-method
#'
#' @description The `MChromatograms` class allows to store
#'     [Chromatogram()] objects in a `matrix`-like
#'     two-dimensional structure.
#'
#' @details The `MChromatograms` class extends the base `matrix` class
#'     and hence allows to store [Chromatogram()] objects in a
#'     two-dimensional array. Each row is supposed to contain
#'     `Chromatogram` objects for one MS data *slice* with a common
#'     m/z and rt range. Columns contain `Chromatogram` objects from the
#'     same sample.
#'
#' @export
#'
#' @return For `[`: the subset of the `MChromatograms` object. If a
#'     single element is extracted (e.g. if `i` and `j` are of length
#'     1) a [Chromatogram()] object is returned. Otherwise (if
#'     `drop = FALSE`, the default, is specified) a `MChromatograms`
#'     object is returned. If `drop = TRUE` is specified, the method
#'     returns a `list` of `Chromatogram` objects.
#'
#'     For `phenoData`: an `AnnotatedDataFrame` representing the
#'     pheno data of the object.
#'
#'     For `pData`: a `data.frame` representing the pheno data of
#'     the object.
#'
#'     For `$`: the value of the corresponding column in the pheno data
#'     table of the object.
#'
#'     For all other methods see function description.
#'
#' @rdname MChromatograms-class
#'
#' @seealso `Chromatogram()] for the class representing chromatogram
#'     data.
#'     [chromatogram()] for the method to extract a `MChromatograms` object
#'     from a `MSnExp` or `OnDiskMSnExp` object.
#'     [readSRMData()` for the function to read chromatographic data
#'     of an SRM/MRM experiment.
#'
#' @author Johannes Rainer
#'
#' @param all for `clean`: `logical(1)` whether all 0-intensities
#'     should be removed (`all = TRUE`), or whether 0-intensities
#'     adjacent to peaks should be kept (`all = FALSE`; default).
#'
#' @param binSize for `bin`: `numeric(1)` with the size of the bins
#'     (in seconds).
#'
#' @param breaks For `bin`: `numeric` defining the bins. Usually not
#'     required as the function calculates the bins automatically based on
#'     `binSize` and the retention time range of chromatograms in the same
#'     row.
#'
#' @param col for `plot`: the color to be used for plotting. Either a
#'     vector of length 1 or equal to `ncol(x)`.
#'
#' @param data for `MChromatograms`: a `list` of [Chromatogram()] objects.
#'
#' @param drop for `[`: `logical(1)` whether to drop the dimensionality of the
#'     returned object (if possible). The default is `drop = FALSE`, i.e. each
#'     subsetting returns a `MChromatograms` object (or a `Chromatogram` object
#'     if a single element is extracted).
#'
#' @param featureData for `MChromatograms`: either a `data.frame` or
#'   `AnnotatedDataFrame` with additional information for each row of
#'   chromatograms.
#'
#' @param fun for `bin`: function to be used to aggregate the intensity
#'     values falling within each bin.
#'
#' @param FUN for `transformIntensity`: function to transform chromatograms'
#'     intensity values. Defaults to `FUN = identity`.
#'
#' @param i for `[`: `numeric`, `logical` or `character`
#'     defining which row(s) to extract.
#'
#' @param intensity for `filterIntensity`: `numeric(1)` or `function` to use to
#'   filter intensities. See description for details.
#'
#' @param j for `[`: `numeric`, `logical` or `character`
#'     defining which columns(s) to extract.
#'
#' @param lty for `plot`: the line type (see `plot` in the `graphics` package
#'     for more details). Can be either a vector of length 1 or of length equal
#'     to `ncol(x)`.
#'
#' @param method `character(1)`. For `normalise`: defining whether each
#'     chromatogram should be normalized to its maximum signal
#'     (`method = "max"`) or total signal (`method = "sum"`).
#'     For `alignRt`: alignment methods (see documentation for `alignRt` in the
#'     [Chromatogram()] help page. Defaults to `method = "closest"`.
#'
#' @param name for `$`, the name of the pheno data column.
#'
#' @param na.rm for `clean`: `logical(1)` whether all `NA`
#'     intensities should be removed prior to clean 0-intensity data points.
#'
#' @param object a `MChromatograms` object.
#'
#' @param phenoData for `MChromatograms`: either a `data.frame`,
#'   `AnnotatedDataFrame` describing the phenotypical information of the
#'   samples.
#'
#' @param type for `plot`: the type of plot (see `plot` from the `graphics`
#'     package for more details). Can be either a vector of length 1 or of
#'     length equal to `ncol(x)`.
#'
#' @param value for `[<-`: the replacement object(s). Can be a `list`
#'     of [Chromatogram()` objects or, if length of `i` and `j` are 1, a
#'     single `Chromatogram` object.
#'
#'     For `pData<-`: a `data.frame` with the number of rows matching
#'     the number of columns of `object`.
#'
#'     For `colnames`: a `character` with the new column names.
#'
#' @param x for all methods: a `MChromatograms` object.
#'
#' @param y for `alignRt`: a [Chromatogram()] object against which `x` should be
#'     aligned against.
#'
#' @param ... for `MChromatograms`: additional parameters to be passed to the
#'     `matrix` constructor, such as `nrow`, `ncol` and `byrow`.
#'     For `compareChromatograms`: ignored.
#'
#' @inheritParams Chromatogram-class
#'
#' @section Object creation:
#'
#' `MChromatograms` are returned by a [chromatogram()] function from an `MSnExp`
#' or `OnDiskMSnExp`. Alternatively, the `MChromatograms` constructor function
#' can be used.
#'
#' @section Data access:
#'
#' - `$` and `$<-`: get or replace individual columns of the object's phenodata.
#'
#' - `colnames` and `colnames<-`: replace or set the column names of the
#'   `MChromatograms` object. Does also set the `rownames` of the `phenoData`.
#'
#' - `fData`: return the feature data as a `data.frame`.
#'
#' - `fData<-`: replace the object's feature data by passing a `data.frame`.
#'
#' - `featureData`: return the feature data.
#'
#' - `featureData<-`: replace the object's feature data.
#'
#' - `featureNames`: returns the feature names of the `MChromatograms` object.
#'
#' - `featureNames<-`: set the feature names.
#'
#' - `fvarLabels`: return the feature data variable names (i.e. column names).
#'
#' - `isEmpty`: returns `TRUE` if the `MChromatograms` object or all of its
#'   `Chromatogram` objects is/are empty or contain only `NA` intensities.
#'
#' - `mz`: returns the m/z for each row of the `MChromatograms` object
#'   as a two-column `matrix` (with columns `"mzmin"` and `"mzmax"`).
#'
#' - `pData`: accesses the phenotypical description of the samples. Returns a
#'   `data.frame`.
#'
#' - `pData<-`: replace the phenotype data.
#'
#' - `phenoData`: accesses the phenotypical description of the samples. Returns
#'   an `AnnotatedDataFrame` object.
#'
#' - `polarity`: returns the polarity of the scans/chromatograms: `1`, `0` or
#'   `-1` for positive, negative or unknown polarity.
#'
#' - `precursorMz`: return the precursor m/z from the chromatograms. The
#'   method returns a `matrix` with 2 columns (`"mzmin"` and `"mzmax"`) and as
#'   many rows as there are rows in the `MChromatograms` object. Each row
#'   contains the precursor m/z of the chromatograms in that row. An error is
#'   thrown if the chromatograms within one row have different precursor m/z
#'   values.
#'
#' - `productMz`: return the product m/z from the chromatograms. The method
#'   returns a `matrix` with 2 columns (`"mzmin"` and `"mzmax"`) and as many
#'   rows as there are rows in the `MChromatograms` object. Each row contains
#'   the product m/z of the chromatograms in that row. An error is thrown if
#'   the chromatograms within one row have different product m/z values.
#'
#' - `rownames<-`: replace the rownames (and featureNames) of the object.
#'
#'
#' @section Data subsetting, combining and filtering:
#'
#' - `[` subset (similar to a `matrix`) by row and column (with parameters `i`
#'   and `j`).
#'
#' - `[<-` replace individual or multiple elements. `value` has to be either a
#'   single `Chromatogram` obhect or a `list` of `Chromatogram` objects.
#'
#' - `c` concatenate (row-wise) `MChromatogram` objects with the
#'   **same number of samples (columns)**.
#'
#' - `filterIntensity`: filter each [Chromatogram()] object within the
#'   `MChromatograms` removing data points with intensities below the user
#'   provided threshold. If `intensity` is a `numeric` value, the returned
#'   chromatogram will only contain data points with intensities > `intensity`.
#'   In addition it is possible to provide a function to perform the filtering.
#'   This function is expected to take the input `Chromatogram` (`object`) and
#'   to return a logical vector with the same length then there are data points
#'   in `object` with `TRUE` for data points that should be kept and `FALSE`
#'   for data points that should be removed. See the `filterIntensity`
#'   documentation in the [Chromatogram()] help page for details and examples.
#'
#'
#' @section Data processing and manipulation:
#'
#' - `alignRt`: align all chromatograms in an `MChromatograms` object against
#'   the chromatogram specified with `y`. See documentation on `alignRt` in the
#'   [Chromatogram()] help page.
#'
#' - `bin`: aggregates intensity values of chromatograms in discrete bins
#'   along the retention time axis. By default, individual `Chromatogram`
#'   objects of one row are binned into the same bins. The function returns a
#'   `MChromatograms` object with binned chromatograms.
#'
#' - `clean`: removes 0-intensity data points. Either all of them
#'   (with `all = TRUE`) or all except those adjacent to non-zero
#'   intensities (`all = FALSE`; default). See [clean()] documentation for more
#'   details and examples.
#'
#' - `compareChromatograms`: calculates pairwise similarity score between
#'   chromatograms in `x` and `y` and returns a similarity matrix with the
#'   number of rows corresponding to the number of chromatograms in `x` and
#'   the number of columns to the number of chromatograms in `y`.
#'   If `y` is missing, a pairwise comparison
#'   is performed between all chromatograms in `x`. See documentation on
#'   `compareChromatograms` in the [Chromatogram()] help page for details.
#'
#' - `normalize`, `normalise`: *normalises* the intensities of a chromatogram by
#'   dividing them either by the maximum intensity (`method = "max"`) or total
#'   intensity (`method = "sum"`) of the chromatogram.
#'
#' - `transformIntensity`: allows to manipulate the intensity values of all
#'   chromatograms using a user provided function. See below for examples.
#'
#' @section Data visualization:
#'
#' - `plot`: plots a `MChromatograms` object. For each row in the object one
#'   plot is created, i.e. all [Chromatogram()] objects in the same row are
#'   added to the same plot. If `nrow(x) > 1` the plot area is split into
#'   `nrow(x)` sub-plots and the chromatograms of one row are plotted in
#'   each.
#'
#' @md
#'
#' @examples
#' ## Creating some chromatogram objects to put them into a MChromatograms object
#' ints <- abs(rnorm(25, sd = 200))
#' ch1 <- Chromatogram(rtime = 1:length(ints), ints)
#' ints <- abs(rnorm(32, sd = 90))
#' ch2 <- Chromatogram(rtime = 1:length(ints), ints)
#' ints <- abs(rnorm(19, sd = 120))
#' ch3 <- Chromatogram(rtime = 1:length(ints), ints)
#' ints <- abs(rnorm(21, sd = 40))
#' ch4 <- Chromatogram(rtime = 1:length(ints), ints)
#'
#' ## Create a MChromatograms object with 2 rows and 2 columns
#' chrs <- MChromatograms(list(ch1, ch2, ch3, ch4), nrow = 2)
#' chrs
#'
#' ## Extract the first element from the second column. Extracting a single
#' ## element always returns a Chromatogram object.
#' chrs[1, 2]
#'
#' ## Extract the second row. Extracting a row or column (i.e. multiple elements
#' ## returns by default a list of Chromatogram objects.
#' chrs[2, ]
#'
#' ## Extract the second row with drop = FALSE, i.e. return a MChromatograms
#' ## object.
#' chrs[2, , drop = FALSE]
#'
#' ## Replace the first element.
#' chrs[1, 1] <- ch3
#' chrs
#'
#' ## Add a pheno data.
#' pd <- data.frame(name = c("first sample", "second sample"),
#'     idx = 1:2)
#' pData(chrs) <- pd
#'
#' ## Column names correspond to the row names of the pheno data
#' chrs
#'
#' ## Access a column within the pheno data
#' chrs$name
#'
#' ## Access the m/z ratio for each row; this will be NA for the present
#' ## object
#' mz(chrs)
#'
#'
#' ## Data visualization
#'
#' ## Create some random Chromatogram objects
#' ints <- abs(rnorm(123, mean = 200, sd = 32))
#' ch1 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 231)
#' ints <- abs(rnorm(122, mean = 250, sd = 43))
#' ch2 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 231)
#' ints <- abs(rnorm(125, mean = 590, sd = 120))
#' ch3 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 542)
#' ints <- abs(rnorm(124, mean = 1200, sd = 509))
#' ch4 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 542)
#'
#' ## Combine into a 2x2 MChromatograms object
#' chrs <- MChromatograms(list(ch1, ch2, ch3, ch4), byrow = TRUE, ncol = 2)
#'
#' ## Plot the second row
#' plot(chrs[2, , drop = FALSE])
#'
#' ## Plot all chromatograms
#' plot(chrs, col = c("#ff000080", "#00ff0080"))
#'
#' ## log2 transform intensities
#' res <- transformIntensity(chrs, log2)
#' plot(res)
NULL

#' @rdname MChromatograms-class
setMethod("show", "MChromatograms", function(object) {
    nr <- nrow(object)
    nc <- ncol(object)
    cat(class(object), " with ",
        nr, ifelse(nr == 1, " row and ", " rows and "),
        nc, ifelse(nc == 1, " column\n", " columns\n"),
        sep = "")
    sumFun <- function(z) {
        paste0("length: ", length(z[[1]]))
    }
    if (nr > 0 && nc > 0) {
        if (nr <= 4) {
            out <- apply(object, MARGIN = c(1, 2), sumFun)
            rownames(out) <- paste0("[", 1:nrow(out), ",]")
        }
        else {
            out <- rbind(
                apply(object[c(1, 2), , drop = FALSE], MARGIN = c(1, 2), sumFun),
                rep(" ... ", ncol(object)),
                apply(object[nrow(object) - c(1, 0), , drop = FALSE],
                      MARGIN = c(1, 2), sumFun)
            )
            rownames(out) <- c("[1,]", "[2,]", "...",
                               paste0("[", c(nrow(object) - c(1, 0)), ",]"))
        }
        rn <- rownames(out)
        out <- rbind(rep("<Chromatogram>", ncol(out)), out)
        rownames(out) <- c("", rn)
        print(out, quote = FALSE, right = TRUE)
    }
    cat("phenoData with", length(varLabels(object@phenoData)), "variables\n")
    cat("featureData with", length(fvarLabels(object)), "variables\n")
})

setAs("matrix", "MChromatograms", function(from) {
    res <- new("MChromatograms")
    res@.Data <- from
    res@phenoData <- annotatedDataFrameFrom(from, byrow = FALSE)
    res@featureData <- annotatedDataFrameFrom(from, byrow = TRUE)
    res
})

#' @rdname MChromatograms-class
setMethod("[", "MChromatograms",
          function(x, i, j, drop = FALSE) {
              if (missing(i) & missing(j))
                  return(x)
              if (missing(i))
                  i <- seq_len(nrow(x))
              if (missing(j))
                  j <- seq_len(ncol(x))
              if (is.logical(i))
                  i <- which(i)
              if (is.logical(j))
                  j <- which(j)
              ## Return a single element as a Chromatogram
              if (length(i) == 1 & length(j) == 1)
                  return(x@.Data[i, j, drop = TRUE][[1]])
              pd <- x@phenoData
              fd <- x@featureData
              xclass <- class(x)
              ## Multiple elements, return type depends on drop.
              x <- x@.Data[i = i, j = j, drop = drop]
              if (!drop) {
                  x <- as(x, xclass)
                  pd <- pd[j, ]
                  ## Drop levels
                  pData(pd) <- droplevels(pData(pd))
                  x@phenoData <- pd
                  fd <- fd[i, ]
                  pData(fd) <- droplevels(pData(fd))
                  x@featureData <- fd
              }
              if (validObject(x))
                  x
          })

#' @rdname MChromatograms-class
setReplaceMethod("[", "MChromatograms",
                 function(x, i, j, value) {
                     if(missing(i) & missing(j))
                         return(x)
                     if (missing(i))
                         i <- seq_len(nrow(x))
                     if (missing(j))
                         j <- seq_len(ncol(x))
                     if (is.logical(i))
                         i <- which(i)
                     if (is.logical(j))
                         j <- which(j)
                     if (is.character(i))
                         stop("Sub-setting rows by name currently not supported")
                     if (is.character(j))
                         j <- match(j, colnames(x))
                     ## check i and j - currently we support only replacing
                     ## single elements, rows or columns.
                     if (length(i) > 1 & length(j) > 1)
                         stop("Only replacement of single elements, rows or ",
                              "columns supported")
                     ## Check value.
                     if (is(value, "Chromatogram"))
                         value <- list(value)
                     res <- vapply(value, FUN = is, FUN.VALUE = logical(1),
                                   "Chromatogram")
                     if(!all(res))
                         stop("All elements in 'value' have to be of type ",
                              "'Chromatogram'")
                     x@.Data[i, j] <- value
                     if (validObject(x))
                         x
                 })

#' @rdname MChromatograms-class
setMethod("plot", signature = signature("MChromatograms"),
          function(x, col = "#00000060", lty = 1, type = "l",
                   xlab = "retention time", ylab = "intensity",
                   main = NULL, ...){
              if (isEmpty(x)) {
                  ## Show a warning and plot an empty plot (issue #249)
                  warning("All chromatograms empty")
                  plot(3, 3, pch = NA, xlab = xlab, ylab = ylab, main = main)
                  text(3, 3, labels = "Empty MChromatograms", col = "red")
              } else {
                  nr <- nrow(x)
                  nc <- ncol(x)
                  ## Initialize plot if we've got more than one row.
                  if (nr > 1) {
                      par(mfrow = c(round(sqrt(nr)), ceiling(sqrt(nr))))
                  }
                  for (i in seq_len(nr)) {
                      if (nc > 1)
                          .plotChromatogramList(x[i, , drop = TRUE], col = col,
                                                lty = lty, type = type,
                                                xlab = xlab, ylab = ylab,
                                                main = main, ...)
                      else
                          plot(x[i, 1], col = col, lty = lty, type = type,
                               xlab = xlab, ylab = ylab, main = main, ...)
                  }
              }
          })

#' @rdname MChromatograms-class
setMethod("phenoData", "MChromatograms", function(object) object@phenoData)

#' @rdname MChromatograms-class
setMethod("pData", "MChromatograms", function(object) pData(phenoData(object)))

#' @rdname MChromatograms-class
setReplaceMethod("pData", c("MChromatograms", "data.frame"),
                 function(object, value) {
                     pData(object@phenoData) <- value
                     if (validObject(object))
                         object
                 })

#' @rdname MChromatograms-class
setMethod("$", "MChromatograms", function(x, name) {
    ## eval(substitute(pData(x)$NAME_ARG, list(NAME_ARG = name)))
    pData(x)[[name]]
})
#' @rdname MChromatograms-class
setReplaceMethod("$", "MChromatograms", function(x, name, value) {
    pData(x)[[name]] <- value
    if(validObject(x))
        x
})

#' @rdname MChromatograms-class
setReplaceMethod("colnames", "MChromatograms", function(x, value) {
    colnames(x@.Data) <- value
    rownames(pData(x)) <- value
    if (validObject(x))
        x
})

#' @rdname MChromatograms-class
setMethod("sampleNames", "MChromatograms", function(object)
    sampleNames(object@phenoData))

#' @rdname MChromatograms-class
setReplaceMethod("sampleNames", "MChromatograms",
                 function(object, value) {
                     colnames(object) <- value
                     object
                 })

#' @rdname MChromatograms-class
setMethod("isEmpty", "MChromatograms", function(x) {
    (nrow(x) == 0 | all(unlist(lapply(x, isEmpty))))
})

#' @rdname MChromatograms-class
setMethod("featureNames", "MChromatograms", function(object)
          featureNames(featureData(object)))

#' @rdname MChromatograms-class
setReplaceMethod("featureNames", "MChromatograms", function(object, value) {
    if (length(value) != nrow(object))
        stop("length of 'value' has to match the number of rows of the ",
             "'MChromatograms' object")
    rownames(object) <- value
    featureNames(featureData(object)) <- value
    if (validObject(object))
        object
})

#' @rdname MChromatograms-class
setMethod("featureData", "MChromatograms", function(object) object@featureData)

#' @rdname MChromatograms-class
setReplaceMethod("featureData", "MChromatograms", function(object, value) {
    if (is.data.frame(value))
        value <- AnnotatedDataFrame(value)
    if (!is(value, "AnnotatedDataFrame"))
        stop("'value' has to be either a 'data.frame' or an ",
             "'AnnotatedDataFrame'")
    if (nrow(value) != nrow(object))
        stop("nrow of 'value' has to match the number of rows of 'object'")
    object@featureData <- value
    rownames(object) <- rownames(value)
    if (validObject(object))
        object
})

#' @rdname MChromatograms-class
setMethod("fData", "MChromatograms", function(object) pData(object@featureData))

#' @rdname MChromatograms-class
setReplaceMethod("fData", "MChromatograms", function(object, value) {
    if (!is.data.frame(value))
        stop("'value' has to be a 'data.frame'")
    if (nrow(value) != nrow(object))
        stop("nrow of 'value' has to match the number of rows of 'object'")
    pData(object@featureData) <- value
    rownames(object) <- rownames(value)
    if (validObject(object))
        object
})

#' @rdname MChromatograms-class
setMethod("fvarLabels", "MChromatograms", function(object)
    varLabels(featureData(object)))

#' @rdname MChromatograms-class
setReplaceMethod("rownames", "MChromatograms", function(x, value) {
    rownames(x@.Data) <- value
    rownames(x@featureData) <- value
    if (validObject(x))
        x
})

#' @rdname MChromatograms-class
setMethod("precursorMz", "MChromatograms", function(object) {
    .mz_chromatograms(object, mz = "precursorMz")
})

#' @rdname MChromatograms-class
setMethod("productMz", "MChromatograms", function(object) {
    .mz_chromatograms(object, mz = "productMz")
})

#' @rdname MChromatograms-class
setMethod("mz", "MChromatograms", function(object) {
    .mz_chromatograms(object, mz = "mz")
})

#' @rdname MChromatograms-class
setMethod("polarity", "MChromatograms", function(object) {
    if (any(fvarLabels(object) == "polarity"))
        fData(object)$polarity
    else
        rep(-1, nrow(object))
})

#' @rdname MChromatograms-class
setMethod("bin", "MChromatograms", .bin_MChromatograms)

setMethod("addProcessing", "MChromatograms", function(object, FUN, ...) {
    if (missing(FUN))
        return(object)
    object@.Data <- matrix(lapply(object, FUN, ...),
                           nrow = nrow(object), dimnames = dimnames(object))
    validObject(object)
    object
})

#' @rdname MChromatograms-class
setMethod("clean", "MChromatograms", function(object, all = FALSE,
                                              na.rm = FALSE) {
    addProcessing(object, FUN = clean, all = all, na.rm = na.rm)
})

#' @rdname MChromatograms-class
setMethod("normalize", "MChromatograms", function(object, method = c("max",
                                                                     "sum")) {
    method <- match.arg(method)
    addProcessing(object, FUN = .normalize_chromatogram, method = method)
})

#' @rdname MChromatograms-class
setMethod("filterIntensity", "MChromatograms", function(object,
                                                        intensity = 0, ...) {
    addProcessing(object, FUN = .filter_intensity_chromatogram,
                  intensity = intensity, ...)
})

#' @rdname MChromatograms-class
setMethod("alignRt", signature = c(x = "MChromatograms", y = "Chromatogram"),
          function(x, y, method = c("closest", "approx"), ...) {
              addProcessing(x, FUN = .align_chromatogram, y = y,
                            method = method, ...)
          })

#' @rdname MChromatograms-class
setMethod("c", "MChromatograms",
          function(x, ...) .bind_rows_chromatograms(x, ...))

#' @rdname MChromatograms-class
setMethod("compareChromatograms",
          signature = c(x = "MChromatograms", y = "missing"),
          function(x, y, ALIGNFUN = alignRt, ALIGNFUNARGS = list(),
                   FUN = cor, FUNARGS = list(use = "pairwise.complete.obs"),
                   ...) {
              .compare_chromatograms_self(
                  x, ALIGNFUN = alignRt, ALIGNFUNARGS = ALIGNFUNARGS,
                  FUN = cor, FUNARGS = FUNARGS)
          })

#' @rdname MChromatograms-class
setMethod("compareChromatograms",
          signature = c(x = "MChromatograms", y = "MChromatograms"),
          function(x, y, ALIGNFUN = alignRt, ALIGNFUNARGS = list(),
                   FUN = cor, FUNARGS = list(use = "pairwise.complete.obs"),
                   ...) {
              .compare_chromatograms(
                  x, y, ALIGNFUN = alignRt, ALIGNFUNARGS = ALIGNFUNARGS,
                  FUN = cor, FUNARGS = FUNARGS)
          })

.compare_chromatograms <- function(x, y, ALIGNFUN = alignRt,
                                   ALIGNFUNARGS = list(), FUN = cor,
                                   FUNARGS = list(
                                       use = "pairwise.complete.obs")) {
    if (inherits(x, "MChromatograms")) {
        if (ncol(x) > 1)
            stop("Currently only single column MChromatograms are supported.")
        x <- unlist(x)
    }
    if (inherits(y, "MChromatograms")) {
        if (ncol(y) > 1)
            stop("Currently only single column MChromatograms are supported.")
        y <- unlist(y)
    }
    nx <- length(x)
    ny <- length(y)

    m <- matrix(NA_real_, nrow = nx, ncol = ny)
    for (i in seq_len(nx)) {
        for (j in seq_len(ny)) {
            m[i, j] <- compareChromatograms(
                x[[i]], y[[j]], ALIGNFUN = ALIGNFUN,
                ALIGNFUNARGS = ALIGNFUNARGS,
                FUN = FUN, FUNARGS = FUNARGS)
        }
    }
    m
}

.compare_chromatograms_self <- function(x, ALIGNFUN = alignRt,
                                        ALIGNFUNARGS = list(), FUN = cor,
                                        FUNARGS = list(
                                            use = "pairwise.complete.obs")) {
    if (inherits(x, "MChromatograms")) {
        if (ncol(x) > 1)
            stop("Currently only single column MChromatograms are supported.")
        x <- unlist(x)
    }
    nx <- length(x)

    m <- matrix(NA_real_, nrow = nx, ncol = nx)
    for (i in seq_len(nx)) {
        for (j in seq_len(nx)) {
            if (i > j)
                next
            m[j, i] <- m[i, j] <- compareChromatograms(
                           x[[i]], x[[j]], ALIGNFUN = ALIGNFUN,
                           ALIGNFUNARGS = ALIGNFUNARGS,
                           FUN = FUN, FUNARGS = FUNARGS)
        }
    }
    m
}

#' @rdname MChromatograms-class
setMethod("transformIntensity", "MChromatograms", function(object,
                                                           FUN = identity) {
    addProcessing(object, FUN = .chromatogram_transform_intensity, FUN)
})
