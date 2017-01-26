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
           max = div <- apply(exprs(object), 1L, max, na.rm = TRUE),
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


##' This function calculates the column-wise coefficient of variation (CV), i.e.
##' the ration between the standard deviation and the mean, for the features
##' in an \code{"\linkS4class{MSnSet}"}. The CVs are calculated for the groups
##' of features defined by \code{groupBy}. For groups defined by single features,
##' \code{NA} is returned.
##'
##' @title Calculates coeffivient of variation for features
##' @param x An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param groupBy An object of class \code{factor} defining how to summerise the features.
##' @param na.rm A \code{logical} defining whether missing values should be removed.
##' @param norm One of 'none' (default), 'sum', 'max', 'center.mean', 'center.median'
##' 'quantiles' or 'quantiles.robust' defining if and how the data should be normalised
##' prior to CV calculation. See \code{\link{normalise}} for more details.
##' @return A \code{matrix} of dimensions \code{length(levels(groupBy))} by \code{ncol(x)}
##' with the respecive CVs.
##' @author Laurent Gatto <lg390@@cam.ac.uk>,
##' Sebastian Gibb <mail@@sebastiangibb.de>
##' @seealso \code{\link{combineFeatures}}
##' @examples
##' data(msnset)
##' msnset <- msnset[1:4]
##' gb <- factor(rep(1:2, each = 2))
##' featureCV(msnset, gb)
featureCV <- function(x, groupBy, na.rm = TRUE,
                      norm = c("sum", "max", "none",
                        "center.mean", "center.median",
                        "quantiles", "quantiles.robust")) {
  norm <- match.arg(norm)
  if (norm != "none")
    x <- normalise(x, method = norm)

  ans <- utils.applyColumnwiseByGroup(exprs(x), groupBy=groupBy,
                                      FUN=function(y, ...) {
                                        utils.colSd(y, ...)/
                                          colMeans(y, ...)}, na.rm=na.rm)
  colnames(ans) <- paste("CV", colnames(ans), sep = ".")
  ans
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

  ans <- utils.applyColumnwiseByGroup(exprs(x), groupBy=groupBy,
                                      FUN=function(y) {
                                        nrow(y)-colSums(is.na(y))})
  colnames(ans) <- sampleNames(x)
  ans
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
    res <- lapply(x, "[", cmn)
    if (!is.null(nms))
        names(res) <- nms
    return(MSnSetList(x = res,
                      log = list(call = match.call())))
}

##' Select feature variables to be retained.
##'
##' @title Select feature variables of interest
##' @param object An \code{MSnSet}.

##' @param graphics A \code{logical} (default is \code{TRUE})
##'     indicating whether a shiny application should be used if
##'     available. Otherwise, a text menu is used. Ignored if \code{k}
##'     is not missing.
##' @param fcol A \code{numeric}, \code{local} or \code{character} of
##'     valid feature variables to be passed directly.
##' @return Updated \code{MSnSet}, containing only selected feature
##'     variables.
##' @author Laurent Gatto
##' @examples
##' library("pRolocdata")
##' data(hyperLOPIT2015)
##' ## 5 first feature variables
##' x <- selectFeatureData(hyperLOPIT2015, fcol = 1:5)
##' fvarLabels(x)
##' \dontrun{
##' ## select via GUI
##' x <- selectFeatureData(hyperLOPIT2015)
##' fvarLabels(x)
##' }
selectFeatureData <- function(object,
                              graphics = TRUE,
                              fcol) {
    if (missing(fcol)) {
        if (graphics) {
            if (!requireNamespace("shiny"))
                warning("The shiny package is required to use the graphical interface.")
            fcol <- .selectShinyFeatureData(object)
        } else fcol <- .selectTextFeatureData(object)
    }
    fData(object) <- fData(object)[, fcol, drop = FALSE]
    object
}


.selectTextFeatureData <- function(object)
    select.list(fvarLabels(object), multiple=TRUE)


.selectShinyFeatureData <- function(object) {
    sel <- fv <- fvarLabels(object)
    on.exit(return(sel))

    ui <- shiny::fluidPage(
        title = 'Examples of DataTables',
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::checkboxGroupInput('vars', 'Feature variables',
                               as.list(fv), selected = sel)),
            shiny::mainPanel(shiny::dataTableOutput('fd'))))

    server <- function(input, output) {
        output$fd <- shiny::renderDataTable({
            sel <<- input$vars
            fData(object)[, input$vars, drop = FALSE]
        })
    }
    app <- list(ui=ui, server=server)
    shiny::runApp(app)
}
