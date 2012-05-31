
normalise.MSnSet <- function(object,method,...) {
  if (method=="vsn") {
    e <- exprs(vsn2(exprs(object)))
  } else if (method=="quantiles") {
    e <- preprocessCore::normalize.quantiles(exprs(object),...)
  } else if (method=="quantiles.robust") {
    e <- preprocessCore::normalize.quantiles.robust(exprs(object),...)
  } else {
    switch(method,
           max = div <- apply(exprs(object),1,max,na.rm=TRUE), 
           sum = div <- rowSums(exprs(object),na.rm=TRUE))
    e <- exprs(object)/div
  }
  rownames(e) <- rownames(exprs(object))
  colnames(e) <- colnames(exprs(object))
  exprs(object) <- e
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Normalised (",method,"): ",
                                              date(),
                                              sep=""))
  object@processingData@normalised <- TRUE
  if (validObject(object))
    return(object)
}

combineMatrixFeatures <- function(matr,    ## matrix
                                   groupBy, ## factor
                                   fun = c("mean",
                                     "median",
                                     "weighted.mean",
                                     "sum",
                                     "medpolish"),
                                  ...,    ## additional arguments to fun
                                  verbose=TRUE) {
  if (is.character(fun)) {
    ## Using a predefined function
    fun <- match.arg(fun)
    if (fun == "medpolish") {
      summarisedFeatures <- by(matr,
                               groupBy,
                               function(x) {
                                 medpol <- medpolish(x, trace.iter = verbose, ...)
                                 return(medpol$overall+medpol$col)
                               })
    } else if (fun == "weighted.mean") {
      ## Expecting 'w' argument
      args <- list(...)
      if (is.null(args$w))
        stop("Expecting a weight parameter 'w' for 'weigthed mean'.")
      w <- args$w
      if (is.null(colnames(matr)))
        colnames(matr) <- paste("X",1:ncol(matr),sep="")
      summarisedFeatures <- apply(matr,2,
                                  function(x) {
                                    .data <- data.frame(x=x,groupBy,w=w)
                                    ddply(.data,
                                          "groupBy",
                                          summarise,
                                          wmn = weighted.mean(x,w))
                                  })
      summarisedFeatures <- do.call(cbind, as.list(summarisedFeatures))
      rn <- summarisedFeatures[,1]
      summarisedFeatures <- summarisedFeatures[, grep("wmn", colnames(summarisedFeatures))]
      colnames(summarisedFeatures) <- colnames(matr)
      rownames(summarisedFeatures) <- rn
      return(summarisedFeatures)
    } else {
      ## using either 'sum', 'mean', 'median'
      summarisedFeatures <- by(matr,
                               groupBy,
                               function(x) apply(x, 2, eval(parse(text = fun)),...))
    }
  } else {
    ## using user-defined function
    summarisedFeatures <- by(matr,
                             groupBy,
                             function(x) apply(x, 2, fun, ...))
  }
    return(do.call(rbind, as.list(summarisedFeatures)))
}


combineFeatures <- function(object,  ## MSnSet
                            groupBy, ## factor
                            fun = c("mean",
                              "median",
                              "weighted.mean",
                              "sum",
                              "medpolish"),
                            ...,    ## additional arguments to fun
                            verbose=TRUE) {
  if (is.character(fun)) 
    fun <- match.arg(fun)
  n1 <- nrow(object)
  ## !! order of features in matRes is defined by the groupBy factor !!
  matRes <- as.matrix(combineMatrixFeatures(exprs(object), groupBy, fun, ..., verbose = verbose))  
  fdata <- fData(object)[!duplicated(groupBy),]
  fdata <- fdata[order(unique(groupBy)),] ## ordering fdata according to groupBy factor
  rownames(matRes) <- rownames(fdata)
  colnames(matRes) <- sampleNames(object)
  exprs(object) <- matRes
  fData(object) <- fdata
  if (is.character(fun)) {
    msg <- paste("Combined ",n1," features into ",
                 nrow(object)," using ",fun,sep="")
  } else {
    msg <- paste("Combined ",n1," features into ",
                 nrow(object)," using user-defined function",sep="")
  }
  object@qual <- object@qual[0,]
  object@processingData@merged <- TRUE
  if (verbose) {
    message(msg)
    ## message("Dropping spectrum-level 'qual' slot.")
  }
  object@processingData@processing <- c(object@processingData@processing,
                                        paste(msg,": ",
                                              date(),
                                              sep=""))
  if (validObject(object))
    return(object)
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
##' features aggregated. In the former case, sums of feautres might
##' be the result of 0 (if no feature was quantified) to \code{n} 
##' (if all \code{topN}'s \code{n} features were quantified) features, 
##' and one might want to rescale the sums based on the number of 
##' non-NA features effectively summed.
##' 
##' @title Count the number of quantitfied features.
##' @param object An instance of class \code{"\linkS4class{MSnSet}"}.
##' @param fcol The feature variable to consider when counting the
##' number of quantified featues.
##' @return A \code{matrix} of dimensions
##' \code{length(levels(factor(fData(object)[, fcol])))} by \code{ncol(object)}
##' of integers.
##' @author Laurent Gatto
##' @examples
##' data(itraqdata)
##' x <- quantify(itraqdata, reporters = iTRAQ4)
##' n <- 2
##' x <- topN(x, groupBy = fData(x)$ProteinAccession, n)
##' m <- nQuants(x, fcol = "ProteinAccession")
##' y <- combineFeatures(x, groupBy = fData(x)$ProteinAccession, fun = sum)
##' stopifnot(dim(n) == dim(y))
##' head(exprs(y))
##' head(exprs(y) * (n/m))
nQuants <- function(object, fcol) {
  .count <- function(x) {
    m <- rep(nrow(x), ncol(x))
    nna <- apply(x, 2, function(.x) sum(is.na(.x)))
    m - nna
  }
  if (class(object) != "MSnSet")
    stop("'object' must be of class 'MSnSet'.")
  if (missing(fcol))
    stop("'fcol' is required.")
  if (!fcol %in% fvarLabels(object))
    stop("'fcol' not found in fvarLabels(object).")
  res <- by(exprs(object),
            factor(fData(object)[, fcol]),
            .count)
  if (ncol(object) == 1) {
    ans <- as.matrix(res)
  } else {
    ans <- do.call(rbind, res)
  }
  colnames(ans) <- sampleNames(object)
  return(ans)
}

