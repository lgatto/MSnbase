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


combineFeatures <- function(object, # MSnSet
                            groupBy, # factor
                            redundancy.handler=c("ignore","unique.only"), # to be exanded to include additional methods
                            ... # further arguments to combine features
                            )
{
   # wrapper to combineFeatures to handle redundancy in feature to factor mapping
   # e.g. peptide-to-protein redundancy
   if(!is.list(groupBy)){
      result <- .combineFeatures(object, groupBy, ...)
   }else{
      # handling of the redundancy
      if(any(names(groupBy) != rownames(fData(object))))
         stop("names of groupBy list do not match fData of the MSnSet object")
      redundancy.handler <- match.arg(redundancy.handler)
      if(redundancy.handler == "ignore"){
         
         expansion.index <- rep(seq_len(nrow(object)), sapply(groupBy, length))
         new.exprs <- exprs(object)[expansion.index,]
         rownames(new.exprs) <- NULL
         groupBy.idx <- sapply(fData(object), identical, groupBy)
         new.feature.data <- fData(object)[expansion.index,]
         new.feature.data[,groupBy.idx] <- unlist(groupBy)
         rownames(new.feature.data) <- NULL
         new.object <- new("MSnSet", exprs = new.exprs, 
                           featureData = new("AnnotatedDataFrame", 
                                             data = new.feature.data),
                           phenoData=phenoData(object))
         result <- .combineFeatures(new.object, unlist(groupBy), ...)
      }else if(redundancy.handler == "unique.only"){
         idx.unique <- sapply(groupBy, length) < 2
         object <- object[idx.unique,]
         groupBy <- unlist(groupBy[idx.unique])
         result <- .combineFeatures(object, groupBy, ...)
      }else{
         stop("Method \"", redundancy.handler, 
              "\" for handing the redundancy is not implemented!", sep='')
      }
      return(result)
   }
}


.combineFeatures <- function(object,  ## MSnSet
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
##' 'quantiles' or 'quantiles.robust' defining if and how the data should be 
