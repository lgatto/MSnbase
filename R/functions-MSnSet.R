
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
                                   fun=c("mean",
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
                                 medpol <- medpolish(x,trace.iter=verbose,...)
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
                                    ddply(.data,"groupBy",summarise,
                                          wmn=weighted.mean(x,w))
                                  })
      summarisedFeatures <- do.call(cbind,summarisedFeatures)
      rn <- summarisedFeatures[,1]
      summarisedFeatures <- summarisedFeatures[,grep("wmn",colnames(summarisedFeatures))]
      colnames(summarisedFeatures) <- colnames(matr)
      rownames(summarisedFeatures) <- rn
      return(summarisedFeatures)
    } else {
      ## using either 'sum', 'mean', 'median'
      summarisedFeatures <- by(matr,
                               groupBy,
                               function(x) apply(x,2,eval(parse(text=fun)),...))
    }
  } else {
    ## using user-defined function
    summarisedFeatures <- by(matr,
                             groupBy,
                             function(x) apply(x,2,fun,...))
  }
    return(do.call(rbind,summarisedFeatures))
}


combineFeatures <- function(object,  ## MSnSet
                            groupBy, ## factor
                            fun=c("mean",
                              "median",
                              "weighted.mean",
                              "sum",
                              "medpolish"),
                            ...,    ## additional arguments to fun
                            verbose=TRUE) {
  n1 <- nrow(object)
  ## !! order of features in matRes is defined by the groupBy factor !!
  matRes <- as.matrix(combineMatrixFeatures(exprs(object),groupBy,fun,...,verbose=verbose))  
  fdata <- fData(object)[!duplicated(groupBy),]
  fdata <- fdata[order(unique(groupBy)),] ## ordering fdata according to groupBy factor
  rownames(matRes) <- rownames(fdata)
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
    warning("Dropping spectrum-level 'qual' slot.")
  }
  object@processingData@processing <- c(object@processingData@processing,
                                        paste(msg,": ",
                                              date(),
                                              sep=""))
  if (validObject(object))
    return(object)
}

