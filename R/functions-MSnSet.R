
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
    e <- sweep(e, 2L, center, check.margin = FALSE)    
  } else if (method == "center.median") {
    e <- exprs(object)
    center <- apply(e, 2L, median, na.rm = TRUE)
    e <- sweep(e, 2L, center, check.margin = FALSE)       
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
                            cv = TRUE,
                            cv.norm = "sum",
                            verbose = TRUE) {
  if (cv) {
    cv.mat <- featureCV(object, groupBy = groupBy,
                        norm = cv.norm)
  }
  if (is.character(fun)) 
    fun <- match.arg(fun)
  n1 <- nrow(object)
  ## !! order of features in matRes is defined by the groupBy factor !!
  matRes <- as.matrix(combineMatrixFeatures(exprs(object),
                                            groupBy, fun, ...,
                                            verbose = verbose)) 
  fdata <- fData(object)[!duplicated(groupBy),] ## takes the first occurences
  fdata <- fdata[order(unique(groupBy)),] ## ordering fdata according to groupBy factor
  rownames(matRes) <- rownames(fdata)
  colnames(matRes) <- sampleNames(object)
  exprs(object) <- matRes
  if (cv)
    fdata <- cbind(fdata, cv.mat)
  fData(object) <- fdata
  if (is.character(fun)) {
    msg <- paste("Combined ", n1, " features into ",
                 nrow(object) ," using ", fun, sep = "")
  } else {
    msg <- paste("Combined ", n1, " features into ",
                 nrow(object), " using user-defined function",
                 sep = "")
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
  ## update feature names according to the groupBy argument
  ## new in version 1.5.9
  fn <- sort(unique(groupBy))
  featureNames(object) <- fn
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
##' @author Laurent Gatto <lg390@@cam.ac.uk>
##' @seealso \code{\link{combineFeatures}}
##' @examples
##' data(itraqdata)
##' m <- quantify(itraqdata[1:4], reporters = iTRAQ4)
##' gb <- factor(rep(1:2, each = 2))
##' featureCV(m, gb)
featureCV <- function(x, groupBy, na.rm = TRUE,
                      norm = c("sum", "max", "none",
                        "center.mean", "center.median",
                        "quantiles", "quantiles.robust")) {
  groupBy <- as.factor(groupBy)
  norm <- match.arg(norm)
  if (norm != "none")
    x <- normalise(x, method = norm)    
  .sd <- function(x, na.rm = na.rm) {
    if (is.matrix(x) | is.data.frame(x)) {
      ans <- apply(x, 2, sd, na.rm = na.rm)
    } else {
      ans <- rep(NA, length(x))
    }
    return(ans)
  }  
  sds <- by(exprs(x), groupBy, .sd, na.rm)  
  mns <- by(exprs(x), groupBy, colMeans)
  stopifnot(all(names(sds) == names(mns)))
  ans <- t(sapply(seq_along(sds), function(i) sds[[i]]/mns[[i]]))
  if (ncol(x) == 1)
    ans <- t(ans)
  rownames(ans) <- names(sds)
  if (is.null(colnames(ans)))
    colnames(ans) <- seq_len(ncol(ans))
  colnames(ans) <- paste("CV", colnames(ans), sep = ".")
  stopifnot(ncol(ans) == ncol(x))
  stopifnot(nrow(ans) == length(levels(groupBy)))
  return(ans)
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

