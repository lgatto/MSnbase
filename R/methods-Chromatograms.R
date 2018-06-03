#' @rdname Chromatograms-class
#'
#' @param object a \code{Chromatograms} object.
setMethod("show", "Chromatograms", function(object) {
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

setAs("matrix", "Chromatograms", function(from) {
    res <- new("Chromatograms")
    res@.Data <- from
    res@phenoData <- as(annotatedDataFrameFrom(from, byrow = FALSE),
                        "NAnnotatedDataFrame")
    res@featureData <- annotatedDataFrameFrom(from, byrow = TRUE)
    res
})

#' @rdname Chromatograms-class
#'
#' @description \code{Chromatograms} objects can, just like a \code{matrix},
#'     be subsetted using the \code{[} method. Single elements, rows or columns
#'     can be replaced using e.g. \code{x[1, 1] <- value} where \code{value}
#'     has to be a \code{Chromatogram} object or a \code{list} of such objects.
#'
#' @note Subsetting with \code{[} will always return a \code{Chromatograms}
#'     object (with the exception of extracting a single element)
#'     unless \code{drop = TRUE} is specified. This is different from the
#'     default subsetting behaviour of \code{matrix}-like objects.
#' 
#' @param x For all methods: a \code{Chromatograms} object.
#' 
#' @param i For \code{[}: \code{numeric}, \code{logical} or \code{character}
#'     defining which row(s) to extract.
#'
#' @param j For \code{[}: \code{numeric}, \code{logical} or \code{character}
#'     defining which columns(s) to extract.
#'
#' @param drop For \code{[}: \code{logical(1)} whether to drop the
#'     dimensionality of the returned object (if possible). The default is
#'     \code{drop = FALSE}, i.e. each subsetting returns a \code{Chromatograms}
#'     object (or a \code{Chromatogram} object if a single element is
#'     extracted).
#'
#' @return For \code{[}: the subset of the \code{Chromatograms} object. If a
#'     single element is extracted (e.g. if \code{i} and \code{j} are of length
#'     1) a \code{\link{Chromatogram}} object is returned. Otherwise (if
#'     \code{drop = FALSE}, the default, is specified) a \code{Chromatograms}
#'     object is returned. If \code{drop = TRUE} is specified, the method
#'     returns a \code{list} of \code{Chromatogram} objects.
#'
#'     For \code{phenoData}: an \code{NAnnotatedDataFrame} representing the
#'     pheno data of the object.
#'
#'     For \code{pData}: a \code{data.frame} representing the pheno data of
#'     the object.
#'
#'     For \code{$}: the value of the corresponding column in the pheno data
#'     table of the object.
setMethod("[", "Chromatograms",
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
              ## Multiple elements, return type depends on drop.
              x <- x@.Data[i = i, j = j, drop = drop]
              if (!drop) {
                  x <- as(x, "Chromatograms")
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

#' @rdname Chromatograms-class
#'
#' @param value For \code{[<-}: the replacement object(s). Can be a \code{list}
#'     of \code{\link{Chromatogram}} objects or, if length of \code{i} and
#'     \code{j} are 1, a single \code{\link{Chromatogram}} object.
#' 
#'     For \code{pData<-}: a \code{data.frame} with the number of rows matching
#'     the number of columns of \code{object}.
#'
#'     For \code{colnames}: a \code{character} with the new column names.
setReplaceMethod("[", "Chromatograms",
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

#' @rdname Chromatograms-class
#' 
#' @description \code{plot}: plots a \code{Chromatograms} object. For each row
#'     in the object one plot is created, i.e. all \code{\link{Chromatogram}}
#'     objects in the same row are added to the same plot.
#'
#' @details \code{plot}: if \code{nrow(x) > 1} the plot area is split into
#'     \code{nrow(x)} sub-plots and the chromatograms of one row are plotted in
#'     each.
#'
#' @param col For \code{plot}: the color to be used for plotting. Either a
#'     vector of length 1 or equal to \code{ncol(x)}.
#'
#' @param lty For \code{plot}: the line type (see \code{\link[graphics]{plot}}
#'     for more details. Can be either a vector of length 1 or of length equal
#'     to \code{ncol(x)}.
#'
#' @param type For \code{plot}: the type of plot (see
#'     \code{\link[graphics]{plot}} for more details. Can be either a vector
#'     of length 1 or of length equal to \code{ncol(x)}.
#' 
#' @inheritParams Chromatogram-class
#'
#' @examples
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
#' ## Combine into a 2x2 Chromatograms object
#' chrs <- Chromatograms(list(ch1, ch2, ch3, ch4), byrow = TRUE, ncol = 2)
#'
#' ## Plot the second row
#' plot(chrs[2, , drop = FALSE])
#'
#' ## Plot all chromatograms
#' plot(chrs, col = c("#ff000080", "#00ff0080"))
setMethod("plot", signature = signature("Chromatograms"),
          function(x, col = "#00000060", lty = 1, type = "l",
                   xlab = "retention time", ylab = "intensity",
                   main = NULL, ...){
              if (isEmpty(x)) {
                  ## Show a warning and plot an empty plot (issue #249)
                  warning("All chromatograms empty")
                  plot(3, 3, pch = NA, xlab = xlab, ylab = ylab, main = main)
                  text(3, 3, labels = "Empty Chromatograms", col = "red")
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

#' @rdname Chromatograms-class
#'
#' @description \code{phenoData}: accesses the phenotypical desccription of the
#'     samples. Returns an \code{NAnnotatedDataFrame} object.
setMethod("phenoData", "Chromatograms", function(object) object@phenoData)

#' @rdname Chromatograms-class
#'
#' @description \code{pData}: accesses the phenotypical description of the
#'     samples. Returns a \code{data.frame}.
setMethod("pData", "Chromatograms", function(object) pData(phenoData(object)))

#' @rdname Chromatograms-class
#'
#' @description \code{pData<-}: replace the phenotype data.
setReplaceMethod("pData", c("Chromatograms", "data.frame"),
                 function(object, value) {
                     pData(object@phenoData) <- value
                     if (validObject(object))
                         object
                 })

#' @rdname Chromatograms-class
#'
#' @description \code{$} and \code{$<-}: get or replace individual columns of
#'     the object's pheno data.
#'
#' @param name For \code{$}, the name of the pheno data column.
setMethod("$", "Chromatograms", function(x, name) {
    ## eval(substitute(pData(x)$NAME_ARG, list(NAME_ARG = name)))
    pData(x)[[name]]
})
#' @rdname Chromatograms-class
setReplaceMethod("$", "Chromatograms", function(x, name, value) {
    pData(x)[[name]] <- value
    if(validObject(x))
        x
})

#' @rdname Chromatograms-class
#'
#' @description \code{colnames<-}: replace or set the column names of the
#'     \code{Chromatograms} object. Does also set the \code{rownames} of the
#'     \code{phenoData}.
setReplaceMethod("colnames", "Chromatograms", function(x, value) {
    colnames(x@.Data) <- value
    rownames(pData(x)) <- value
    if (validObject(x))
        x
})

#' @rdname Chromatograms-class
#'
#' @description \code{sampleNames}: get the sample names.
setMethod("sampleNames", "Chromatograms", function(object)
    sampleNames(object@phenoData))

#' @rdname Chromatograms-class
#'
#' @description \code{sampleNames<-}: replace or set the sample names of the
#'     \code{Chromatograms} object (i.e. the \code{rownames} of the pheno data
#'     and \code{colnames} of the data matrix.
setReplaceMethod("sampleNames", "Chromatograms",
                 function(object, value) {
                     colnames(object) <- value
                     object
                 })

#' @rdname Chromatograms-class
#'
#' @description \code{isEmpty}: returns \code{TRUE} if the \code{Chromatograms}
#'     object or all of its \code{Chromatogram} objects is/are empty or contain
#'     only \code{NA} intensities.
setMethod("isEmpty", "Chromatograms", function(x) {
    (nrow(x) == 0 | all(unlist(lapply(x, isEmpty))))
})

#' @rdname Chromatograms-class
#' 
#' @description \code{featureNames}: returns the feature names of the
#'     \code{Chromatograms} object.
setMethod("featureNames", "Chromatograms", function(object)
          featureNames(featureData(object)))

#' @rdname Chromatograms-class
#'
#' @description \code{featureNames<-}: set the feature names.
setReplaceMethod("featureNames", "Chromatograms", function(object, value) {
    if (length(value) != nrow(object))
        stop("length of 'value' has to match the number of rows of the ",
             "'Chromatograms' object")
    rownames(object) <- value
    featureNames(featureData(object)) <- value
    if (validObject(object))
        object
})

#' @rdname Chromatograms-class
#'
#' @description \code{featureData}: return the feature data.
setMethod("featureData", "Chromatograms", function(object) object@featureData)

#' @rdname Chromatograms-class
#'
#' @description \code{featureData<-}: replace the object's feature data.
setReplaceMethod("featureData", "Chromatograms", function(object, value) {
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

#' @rdname Chromatograms-class
#'
#' @description \code{fData}: return the feature data as a \code{data.frame}.
setMethod("fData", "Chromatograms", function(object) pData(object@featureData))

#' @rdname Chromatograms-class
#'
#' @description \code{fData<-}: replace the object's feature data by passing a
#'     \code{data.frame}
setReplaceMethod("fData", "Chromatograms", function(object, value) {
    if (!is.data.frame(value))
        stop("'value' has to be a 'data.frame'")
    if (nrow(value) != nrow(object))
        stop("nrow of 'value' has to match the number of rows of 'object'")
    pData(object@featureData) <- value
    rownames(object) <- rownames(value)
    if (validObject(object))
        object
})

#' @rdname Chromatograms-class
#'
#' @description \code{fvarLabels}: return the feature data variable names (i.e.
#'     column names).
setMethod("fvarLabels", "Chromatograms", function(object)
    varLabels(featureData(object)))

#' @rdname Chromatograms-class
#'
#' @description \code{rownames<-}: replace the rownames (and featureNames) of
#'     the object.
setReplaceMethod("rownames", "Chromatograms", function(x, value) {
    rownames(x@.Data) <- value
    rownames(x@featureData) <- value
    if (validObject(x))
        x
})

#' @rdname Chromatograms-class
#'
#' @description
#'
#' \code{precursorMz}: return the precursor m/z from the chromatograms. The
#' method returns a \code{matrix} with 2 columns (\code{"mzmin"} and
#' \code{"mzmax"}) and as many rows as there are rows in the
#' \code{Chromatograms} object. Each row contains the precursor m/z of the
#' chromatograms in that row. An error is thrown if the chromatograms within one
#' row have different precursor m/z values.
setMethod("precursorMz", "Chromatograms", function(object) {
    .mz_chromatograms(object, mz = "precursorMz")
})

#' @rdname Chromatograms-class
#'
#' @description
#'
#' \code{productMz}: return the product m/z from the chromatograms. The
#' method returns a \code{matrix} with 2 columns (\code{"mzmin"} and
#' \code{"mzmax"}) and as many rows as there are rows in the
#' \code{Chromatograms} object. Each row contains the product m/z of the
#' chromatograms in that row. An error is thrown if the chromatograms within one
#' row have different product m/z values.
setMethod("productMz", "Chromatograms", function(object) {
    .mz_chromatograms(object, mz = "productMz")
})

#' @rdname Chromatograms-class
#'
#' @description
#'
#' \code{mz}: returns the m/z for each row of the \code{Chromatograms} object
#' as a two-column \code{matrix} (with columns \code{"mzmin"} and
#' \code{"mzmax"}).
setMethod("mz", "Chromatograms", function(object) {
    .mz_chromatograms(object, mz = "mz")
})

#' @rdname Chromatograms-class
#'
#' @description
#'
#' \code{polarity}: returns the polarity of the scans/chromatograms: `1`,
#' `0` or `-1` for positive, negative or unknown polarity.
setMethod("polarity", "Chromatograms", function(object) {
    if (any(fvarLabels(object) == "polarity"))
        fData(object)$polarity
    else
        rep(-1, nrow(object))
})

#' @description
#'
#' \code{bin} aggregates intensity values of chromatograms in discrete bins
#' along the retention time axis. By default, individual \code{Chromatogram}
#' objects of one row are binned into the same bins. The function returns a
#' \code{Chromatograms} object with binned chromatograms.
#'
#' @param binSize for \code{bin}: \code{numeric(1)} with the size of the bins
#'     (in seconds).
#'
#' @param breaks for \code{bin}: \code{numeric} defining the bins. Usually not
#'     required as the function calculates the bins automatically based on
#'     \code{binSize} and the retention time range of chromatograms in the same
#'     row.
#'
#' @param fun for \code{bin}: function to be used to aggregate the intensity
#'     values falling within each bin.
#'
#' @rdname Chromatograms-class
setMethod("bin", "Chromatograms", .bin_Chromatograms)
