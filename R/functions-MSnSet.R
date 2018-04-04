MSnSet <- function(exprs, fData, pData, ...) {
  if (class(fData) == "data.frame")
    fData <- new("AnnotatedDataFrame", data = fData)
  if (class(pData) == "data.frame")
    pData <- new("AnnotatedDataFrame", data = pData)
  ans <- new("MSnSet",
             exprs = exprs,
             featureData = fData,
             phenoData = pData)
  if (validObject(ans))
      return(ans)
}

normalise_MSnSet <- function(object, method, ...) {
  if (method == "vsn") {
    e <- exprs(vsn2(exprs(object), ...))
  } else if (method == "quantiles") {
    e <- preprocessCore::normalize.quantiles(exprs(object), ...)
  } else if (method == "quantiles.robust") {
    e <- preprocessCore::normalize.quantiles.robust(exprs(object), ...)
  } else if (method == "center.mean") {
    e <- exprs(object)
    center <- colMeans(e, na.rm = TRUE)
    e <- sweep(e, 2L, center, check.margin = FALSE, ...)
  } else if (method == "center.median") {
    e <- exprs(object)
    center <- apply(e, 2L, median, na.rm = TRUE)
    e <- sweep(e, 2L, center, check.margin = FALSE, ...)
  } else if (method == "diff.median") {
      e <- exprs(object)
      med <- median(as.numeric(e), na.rm = TRUE)
      cmeds <- apply(e, 2L, median, na.rm = TRUE)
      e <- sweep(e, 2L, cmeds - med)
  } else {
    switch(method,
           max = div <- .rowMaxs(exprs(object), na.rm = TRUE),
           sum = div <- rowSums(exprs(object), na.rm = TRUE))
    e <- exprs(object)/div
  }
  rownames(e) <- rownames(exprs(object))
  colnames(e) <- colnames(exprs(object))
  exprs(object) <- e
  object@processingData@processing <-
    c(object@processingData@processing,
      paste("Normalised (", method ,"): ",
            date(), sep = ""))
  object@processingData@normalised <- TRUE
  if (validObject(object))
    return(object)
}


##' This function calculates the column-wise coefficient of variation
##' (CV), i.e.  the ration between the standard deviation and the
##' mean, for the features in an [`MSnSet`]. The CVs are calculated
##' for the groups of features defined by `groupBy`. For groups
##' defined by single features, `NA` is returned.
##'
##' @title Calculates coeffivient of variation for features
##' @param x An instance of class [`MSnSet`].
##' @param groupBy An object of class `factor` defining how to
##'     summarise the features.
##' @param na.rm A `logical(1)` defining whether missing values should
##'     be removed.
##' @param norm One of normalisation methods applied prior to CV
##'     calculation. See [normalise()] for more details. Here, the
##'     default is `'none'`, i.e. no normalisation.
##' @param suffix A `character(1)` to be used to name the new CV
##'     columns. Default is `NULL` to ignore this. This argument
##'     should be set when CV values are already present in the
##'     [`MSnSet`] feature variables.
##' @return A `matrix` of dimensions `length(levels(groupBy))` by
##'     `ncol(x)` with the respecive CVs. The column names are formed
##'     by pasting `CV.` and the sample names of object `x`, possibly
##'     suffixed by `.suffix`.
##' @author Laurent Gatto and Sebastian Gibb
##' @seealso [combineFeatures()]
##' @md
##' @examples
##' data(msnset)
##' msnset <- msnset[1:4]
##' gb <- factor(rep(1:2, each = 2))
##' featureCV(msnset, gb)
##' featureCV(msnset, gb, suffix = "2")
featureCV <- function(x, groupBy, na.rm = TRUE,
                      norm = "none",
                      suffix = NULL) {
  if (norm != "none")
    x <- normalise(x, method = norm)

  cv <- rowsd(exprs(x), group = groupBy, reorder = TRUE, na.rm = na.rm) /
      rowmean(exprs(x), group = groupBy, reorder = TRUE, na.rm = na.rm)
  colnames(cv) <- paste("CV", colnames(cv), sep = ".")
  if (!is.null(suffix))
      colnames(cv) <- paste(colnames(cv), suffix, sep = ".")
  cv
}

updateFvarLabels <- function(object, label, sep = ".") {
  if(missing(label))
    label <- getVariableName(match.call(), "object")
  fvarLabels(object) <- paste(fvarLabels(object),
                              label, sep = sep)
  object
}


updateSampleNames <- function(object, label, sep = ".") {
  if (missing(label))
    label <- getVariableName(match.call(), "object")
  sampleNames(object) <- paste(sampleNames(object),
                               label, sep = sep)
  object
}

updateFeatureNames <- function(object, label, sep = ".") {
  if (missing(label))
    label <- getVariableName(match.call(), "object")
  featureNames(object) <- paste(featureNames(object),
                                label, sep = sep)
  object
}

##' This function counts the number of quantified features, i.e
##' non NA quantitation values, for each group of features
##' for all the samples in an \code{"\linkS4class{MSnSet}"} object.
##' The group of features are defined by a feature variable names, i.e
##' the name of a column of \code{fData(object)}.
##'
##' This function is typically used after \code{\link{topN}} and before
##' \code{\link{combineFeatures}}, when the summerising function is
##' \code{sum}, or any function that does not normalise to the number of
##' features aggregated. In the former case, sums of features might
##' be the result of 0 (if no feature was quantified) to \code{n}
##' (if all \code{topN}'s \code{n} features were quantified) features,
##' and one might want to rescale the sums based on the number of
##' non-NA features effectively summed.
##'
##' @title Count the number of quantitfied features.
##' @param x An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param groupBy An object of class \code{factor} defining how to
##'     summerise the features. (Note that this parameter was
##'     previously named \code{fcol} and referred to a feature
##'     variable label. This has been updated in version 1.19.12 for
##'     consistency with other functions.)
##' @return A \code{matrix} of dimensions
##'     \code{length(levels(groupBy))} by \code{ncol(x)}
##' @return A \code{matrix} of dimensions
##'     \code{length(levels(factor(fData(object)[, fcol])))} by
##'     \code{ncol(object)} of integers.
##' @author Laurent Gatto <lg390@@cam.ac.uk>, Sebastian Gibb
##'     <mail@@sebastiangibb.de>
##' @examples
##' data(msnset)
##' n <- 2
##' msnset <- topN(msnset, groupBy = fData(msnset)$ProteinAccession, n)
##' m <- nQuants(msnset, groupBy = fData(msnset)$ProteinAccession)
##' msnset2 <- combineFeatures(msnset,
##'                            groupBy = fData(msnset)$ProteinAccession,
##'                            fun = sum)
##' stopifnot(dim(n) == dim(msnset2))
##' head(exprs(msnset2))
##' head(exprs(msnset2) * (n/m))
nQuants <- function(x, groupBy) {
  if (class(x) != "MSnSet")
    stop("'x' must be of class 'MSnSet'.")

  e <- !is.na(exprs(x))
  mode(e) <- "numeric"
  rowsum(e, group=groupBy, reorder=TRUE)
}

##' Subsets \code{MSnSet} instances to their common feature names.
##'
##' @title Keep only common feature names
##' @param x An instance of class \code{\linkS4class{MSnSet}} or a
##'     \code{list} or \code{MSnSetList} with at least 2 \code{MSnSet}
##'     objects.
##' @param y An instance of class \code{\linkS4class{MSnSet}}. Ignored
##'     if \code{x} is a \code{list}/\code{MSnSetList}.
##' @return An \code{linkS4class{MSnSetList}} composed of the input
##'     \code{MSnSet} containing only common features in the same
##'     order. The names of the output are either the names of the
##'     \code{x} and \code{y} input variables or the names of \code{x}
##'     if a list is provided.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(tan2009r1)
##' data(tan2009r2)
##' cmn <- commonFeatureNames(tan2009r1, tan2009r2)
##' names(cmn)
##' ## as a named list
##' names(commonFeatureNames(list(a = tan2009r1, b = tan2009r2)))
##' ## without message
##' suppressMessages(cmn <- commonFeatureNames(tan2009r1, tan2009r2))
##' ## more than 2 instance
##' data(tan2009r3)
##' cmn <- commonFeatureNames(list(tan2009r1, tan2009r2, tan2009r3))
##' length(cmn)
commonFeatureNames <- function(x, y) {
    if (inherits(x, "MSnSetList"))
        x <- msnsets(x)
    if (inherits(x, "MSnSet")) {
        stopifnot(inherits(y, "MSnSet"))
        nms <- c(getVariableName(match.call(), "x"),
                 getVariableName(match.call(), "y"))
        x <- list(x, y)
        names(x) <- nms
    }
    nms <- names(x)
    cmn <- Reduce(intersect, lapply(x, featureNames))
    message(paste(length(cmn), "features in common"))
    res <- lapply(x, function(xx) xx[cmn, ])
    if (!is.null(nms))
        names(res) <- nms
    return(MSnSetList(x = res,
                      log = list(call = match.call())))
}

##' This function computes the number of features in the group defined
##' by the feature variable \code{fcol} and appends this information
##' in the feature data of \code{object}.
##'
##' @title How many features in a group?
##' @param object An instance of class \code{MSnSet}.
##' @param fcol Feature variable defining the feature grouping
##'     structure.
##' @return An updated \code{MSnSet} with a new feature variable
##'     \code{fcol.nFeatures}.
##' @author Laurent Gatto
##' @examples
##' library(pRolocdata)
##' data("hyperLOPIT2015ms3r1psm")
##' hyperLOPIT2015ms3r1psm <- nFeatures(hyperLOPIT2015ms3r1psm,
##'                                     "Protein.Group.Accessions")
##' i <- c("Protein.Group.Accessions", "Protein.Group.Accessions.nFeatures")
##' fData(hyperLOPIT2015ms3r1psm)[1:10, i]
nFeatures <- function(object, fcol) {
    stopifnot(inherits(object, "MSnSet"))
    stopifnot(fcol %in% fvarLabels(object))
    fcol2 <- paste0(fcol, ".nFeatures")
    if (fcol2 %in% fvarLabels(object))
        stop("'", fcol2, "' already present.")
    k <- table(fData(object)[, fcol])
    k <- k[as.character(fData(object)[, fcol])]
    fData(object)[, fcol2] <- k
    if (validObject(object))
        return(object)
}
