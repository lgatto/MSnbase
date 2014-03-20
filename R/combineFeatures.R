combineFeatures <- function(object, 
                            groupBy, 
                            fun = c("mean",
                                "median",
                                "weighted.mean",
                                "sum",
                                "medpolish"),
                            redundancy.handler = c("unique", "multiple"),
                            cv = TRUE,
                            cv.norm = "sum",
                            verbose = TRUE,
                            ... ## further arguments to fun
                            ) {
    if (is.character(fun))
        fun <- match.arg(fun)
    if (is.list(groupBy)) {
        if (length(groupBy) != nrow(object))
            stop("'length(groupBy)' must be equal to 'nrow(object)': ",
                 length(groupBy), " != ", nrow(object), ".")
        if (!is.null(names(groupBy))) {          
            if (!all(names(groupBy) %in% featureNames(object)))
                stop("'groupBy' names and 'featureNames(object)' don't match.")
            if (!all(names(groupBy) == featureNames(object))) {
                warning("Re-ordering groupBy to match feature names.")
                groupBy <- groupBy[featureNames(object)]
            }
        }
        redundancy.handler <- match.arg(redundancy.handler)         
        result <- combineFeaturesL(object, groupBy, fun,
                                   redundancy.handler,
                                   cv, cv.norm, verbose, ...)        
    } else { ## factor, numeric or character
        result <- combineFeaturesV(object, groupBy, fun,
                                   cv, cv.norm, verbose, ...)
    }
    if (validObject(result))
        return(result)
}

combineFeaturesL <- function(object,   ## MSnSet
                             groupBy,  ## list
                             fun,
                             redundancy.handler,
                             cv = TRUE,
                             cv.norm = "sum",
                             verbose = TRUE,
                             ...    ## additional arguments to fun
                             ) {
    ## handling of the redundancy
    if (redundancy.handler == "multiple") { 
        expansion.index <- rep(seq_len(nrow(object)), sapply(groupBy, length))
        new.exprs <- exprs(object)[expansion.index, ]
        rownames(new.exprs) <- NULL
        new.feature.data <- fData(object)[expansion.index, , drop = FALSE]
        rownames(new.feature.data) <- NULL
        ## TODO: check this
        object <- new("MSnSet", exprs = new.exprs, 
                      featureData = new(
                          "AnnotatedDataFrame", 
                          data = new.feature.data),
                      phenoData = phenoData(object))
        groupBy <- unlist(groupBy)
    } else { ## redundancy.handler == "unique" - checked by match.arg
        idx.unique <- sapply(groupBy, length) < 2
        object <- object[idx.unique, ]
        groupBy <- unlist(groupBy[idx.unique])
    }
    combineFeaturesV(object, groupBy, fun, cv, cv.norm, verbose, ...)
}


combineFeaturesV <- function(object,   ## MSnSet
                             groupBy,  ## factor, character or numeric
                             fun,
                             cv = TRUE,
                             cv.norm = "sum",
                             verbose = TRUE,
                             ...    ## additional arguments to fun
                             ) {
    groupBy <- as.character(groupBy)
    if (cv) {
        cv.mat <- featureCV(object, groupBy = groupBy,
                            norm = cv.norm)
    }
    n1 <- nrow(object)
    ## !! order of features in matRes is defined by the groupBy factor !!
    matRes <- as.matrix(combineMatrixFeatures(exprs(object),
                                              groupBy, fun, 
                                              verbose = verbose,
                                              ...)) 
    fdata <- fData(object)[!duplicated(groupBy), , drop = FALSE] ## takes the first occurences
    fdata <- fdata[order(unique(groupBy)), , drop = FALSE] ## ordering fdata according to groupBy factor
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

combineMatrixFeatures <- function(matr,    ## matrix
                                  groupBy, ## char/factor
                                  fun,
                                  verbose=TRUE, ## additional arguments to fun
                                  ...) {
  if (is.character(fun)) {
    ## Using a predefined function
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
        colnames(matr) <- paste0("X", 1:ncol(matr))
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
                               function(x) apply(x, 2, eval(parse(text = fun)), ...))
    }
  } else {
    ## using user-defined function
    summarisedFeatures <- by(matr,
                             groupBy,
                             function(x) apply(x, 2, fun, ...))
  }
    return(do.call(rbind, as.list(summarisedFeatures)))
}
