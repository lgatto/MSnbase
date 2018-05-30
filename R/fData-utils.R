##' The `expandFeatureVars` and `mergeFeatureVars` respectively expand
##' and merge groups of feature variables. Using these functions, a
##' set of columns in a feature data can be merged into a single new
##' data.frame-column variables and a data.frame-column can be
##' expanded into single feature columns. The original feature
##' variables are removed.
##'
##' @title Expand or merge feature variables
##' @param x An object of class `MSnSet`.
##' @param fcol A `character()` of feature variables to expand (for
##'     `expandFeatureVars`) or merge (for `mergeFeatureVars`).
##' @param prefix A `character(1)` to use as prefix to the new feature
##'     variables. If missing (default), then `fcol` is used
##'     instead. If `NULL`, then no prefix is used.
##' @return An `MSnSet` for expanded (merged) feature variables.
##' @author Laurent Gatto
##' @md
##' @rdname fData-utils
##' @examples
##' library("pRolocdata")
##' data(hyperLOPIT2015)
##' fvarLabels(hyperLOPIT2015)
##' ## Let's merge all svm prediction feature variables
##' (k <- grep("^svm", fvarLabels(hyperLOPIT2015), value = TRUE))
##' hl <- mergeFeatureVars(hyperLOPIT2015, fcol = k, fcol2 = "SVM")
##' fvarLabels(hl)
##' head(fData(hl)$SVM)
##'
##' ## Let's expand the new SVM into individual columns
##' hl2 <- expandFeatureVars(hl, "SVM")
##' fvarLabels(hl2)
##' ## We can set the prefix manually
##' hl2 <- expandFeatureVars(hl, "SVM", prefix = "Expanded")
##' fvarLabels(hl2)
##' ## If we don't want any prefix
##' hl2 <- expandFeatureVars(hl, "SVM", prefix = NULL)
##' fvarLabels(hl2)
expandFeatureVars <- function(x, fcol, prefix) {
    stopifnot(is.character(fcol))
    stopifnot(is.data.frame(fData(x)))
    stopifnot(fcol %in% fvarLabels(x))
    stopifnot(as.logical(ncol(fd <- fData(x)[, fcol])))
    if (missing(prefix)) prefix <- fcol
    for (f in names(fd)) {
        if (is.null(prefix)) f2 <- f
        else f2 <- paste(prefix, f, sep = ".")
        if (f2 %in% fvarLabels(x))
            stop("Feature variable ", f2, " already exists.")
        fData(x)[, f2] <- fd[, f]
    }
    fData(x)[, fcol] <- NULL
    if (validObject(x))
        x
}


##' @param fcol2 A `character(1)` defining the name of the new feature
##'     variable.
##' @md
##' @rdname fData-utils
mergeFeatureVars <- function(x, fcol, fcol2) {
    stopifnot(is.character(fcol2))
    if (fcol2 %in% fvarLabels(x))
        stop("Feature variable ", fcol2, "already exists.")
    stopifnot(is.data.frame(fData(x)))
    if (is.numeric(fcol))
        fcol <- fvarLabels(x)[fcol]
    stopifnot(is.character(fcol) & !anyNA(fcol))
    stopifnot(all(fcol %in% fvarLabels(x)))
    fd <- fData(x)[, fcol]
    fData(x)[[fcol2]] <- fd
    for (f in fcol)
        fData(x)[, f] <- NULL
    if (validObject(x))
        x
}
