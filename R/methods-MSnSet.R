#######################################################################
## NOTE: MSnSet is based in the ExpressionSet class defined in the
##       Biobase package and should have the exact same functionnality.
##       Some of the code is heavily inspired and sometimes
##       directly copies from the respective ExpressionSet method.

setMethod("initialize", "MSnSet",
          function(.Object,
                   assayData,
                   phenoData,
                   featureData,
                   experimentData,
                   exprs = new("matrix"),
                   ... ) {
            if (missing(assayData)) {
              if (missing(phenoData))
                phenoData <- annotatedDataFrameFrom(exprs, byrow=FALSE)
              if (missing(featureData))
                featureData <- annotatedDataFrameFrom(exprs, byrow=TRUE)
              if (missing(experimentData))
                experimentData <- new("MIAPE")
              .Object <- callNextMethod(.Object,
                                        phenoData = phenoData,
                                        featureData = featureData,
                                        exprs = exprs,
                                        experimentData = experimentData,
                                        ...)
            } else if (missing(exprs)) {
              if (missing(phenoData))
                phenoData <- annotatedDataFrameFrom(assayData, byrow=FALSE)
              if (missing(featureData))
                featureData <- annotatedDataFrameFrom(assayData, byrow=TRUE)
              if (missing(experimentData))
                experimentData <- new("MIAPE")
              .Object <- callNextMethod(.Object,
                                        assayData = assayData,
                                        phenoData = phenoData,
                                        featureData = featureData,
                                        experimentData = experimentData,
                                        ...)
            } else stop("provide at most one of 'assayData' or 'exprs' to initialize MSnSet",
                        call.=FALSE)
            Biobase:::.harmonizeDimnames(.Object)
          })



setValidity("MSnSet", function(object) {
  msg <- validMsg(NULL, Biobase:::isValidVersion(object, "MSnSet"))
  msg <- validMsg(msg, Biobase::assayDataValidMembers(assayData(object), c("exprs")))
  if ( nrow(qual(object)) != 0 ) {
    nrow.obs <- nrow(qual(object))
    nrow.exp <- nrow(object)*length(object$reporters)
    if (nrow.obs != nrow.exp)
      msg <- validMsg(msg,
                      "number of rows in assayData and qual slots do not match.")
  }
  if (!inherits(experimentData(object),"MIAPE"))
    msg <- validMsg(msg, 
                    "experimentData slot in MSnSet must be 'MIAPE' object")
  if (is.null(msg)) TRUE else msg
})

setMethod("exprs", signature(object="MSnSet"),
          function(object) assayDataElement(object,"exprs"))

setReplaceMethod("exprs", signature(object="MSnSet",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "exprs", value))

setMethod("show","MSnSet",
          function(object) {
            callNextMethod()
            show(processingData(object))
            invisible(NULL)
          })


setMethod("normalise","MSnSet",
          function(object,method=c("sum","max",
                            "quantiles",
                            "quantiles.robust",
                            "vsn"),...)
          normalise.MSnSet(object,match.arg(method),...)
          )

setMethod("normalize","MSnSet",
          function(object,method,...) normalise(object,method,...)
          )

setMethod("purityCorrect",
          signature=signature("MSnSet","matrix"),
          function(object,impurities) {
            if (ncol(impurities)!=nrow(impurities))
              stop("Impurity matrix must be a square matrix")
            if (ncol(object)!=ncol(impurities))
              stop("Impurity matrix should be ",ncol(object)," by ",ncol(object))
            .purcor <- function(x,.impurities=impurities) {
              keep <- !is.na(x)
              if (sum(keep)>1) 
                x[keep] <- solve(.impurities[keep,keep],x[keep])
              x[x<0] <- NA
              return(x)
            }
            corr.exprs <- apply(exprs(object),1,.purcor)
            exprs(object) <- t(corr.exprs)
            object@processingData@processing <-             
              c(object@processingData@processing,
                paste("Purity corrected: ",date(),sep=""))
            if (validObject(object))
              return(object)
          })


setMethod("dim","MSnSet",function(x) dim(exprs(x)))
setMethod("qual","MSnSet", function(object) object@qual)
## Not sure about these...
## setReplaceMethod("featureNames",
##                  signature(object="MSnSet",
##                            value="character"),
##                  function(object, value) {
##                    object@features = value
##                    if (validObject(object))
##                      return(object)
##                  })

## No proteomicsData anymore (since version 0.2.0 of MSnSet and MSnbase).
## experimentData is not proper MIAPE
## setMethod("proteomicsData","MSnSet",function(object) object@proteomicsData)
## setReplaceMethod("proteomicsData",
##                  signature(object="MSnSet",
##                            value="MIAPE"),
##                  function(object, value) {
##                    object@proteomicsData = value
##                    if (validObject(object))
##                      return(object)
##                  })

setMethod("fileNames",
          signature(object="MSnSet"),
          function(object) processingData(object)@files)

setMethod("processingData",
          signature(object="MSnSet"),
          function(object) object@processingData)

setMethod("msInfo","MSnSet",
          function(object) msInfo(experimentData(object)))

setMethod("meanSdPlot",
          signature="MSnSet",
          definition =
          function(x, ranks=TRUE, xlab = ifelse(ranks, "rank(mean)", "mean"),
                   ylab = "sd", pch  = ".", plot = TRUE, ...)
          vsn::meanSdPlot(exprs(x), ranks=ranks, xlab=xlab, ylab=ylab, pch=pch, plot=plot, ...))

t.MSnSet <- function(x) {
  x@processingData@processing <-             
    c(x@processingData@processing,
      paste("MSnSet transposed: ",date(),sep=""))
  message("Dropping protocolData.")  
  return(new("MSnSet",
             exprs = t(exprs(x)),
             phenoData = featureData(x),
             featureData = phenoData(x),
             experimentData = experimentData(x),
             processingData = processingData(x),
             annotation = annotation(x)))
}


setMethod("[", "MSnSet", function(x, i, j, ...) {
  .Object <- callNextMethod(...)
  ## subsetting qual - requires pData(x)$mz!
  fn <- featureNames(.Object)
  reps <- match(.Object$mz,x$mz)
  qrows <- paste(rep(fn,each=length(reps)),reps,sep=".")
  .Object@qual <- .Object@qual[qrows,]
  .Object@processingData@processing <- c(.Object@processingData@processing,
                                         ifelse(missing(j),
                                                paste("Features subsetted: ",date(),sep=""),
                                                paste("Samples subsetted: ",date(),sep="")))
  if (validObject(.Object))
    return(.Object)
})

  
setAs("MSnSet", "ExpressionSet",
      function (from)
      new("ExpressionSet",
          exprs=exprs(from),
          phenoData=phenoData(from),
          featureData=featureData(from),
          annotation=annotation(from),
          protocolData=protocolData(from))
      )

as.ExpressionSet.MSnSet <- function(x) as(x,"ExpressionSet")

setMethod("write.exprs",
          signature(x="MSnSet"),
          function(x,
                   fDataCols=NULL,
                   file="tmp.txt", quote=FALSE,
                   sep="\t", col.names=NA, ...) {
            res <- exprs(x)
            if (!is.null(fDataCols))
              res <- cbind(res,fData(x)[,fDataCols])
            write.table(res, file=file, quote=quote, sep=sep, col.names=col.names, ...)
          })

setReplaceMethod("experimentData",
                 signature = signature(
                   object = "MSnSet",
                   value = "MIAPE"),
                 function(object, value) {
                   if (!validObject(value))
                     stop("Not a valid MIAPE instance.")
                   object@experimentData <- value
                   object
                 })



setMethod("combine",
          signature=signature(
            x="MSnSet", y="MSnSet"),
          function(x, y, ...) {
            if (class(x) != class(y))
              stop(paste("objects must be the same class, but are ",
                         class(x), ", ", class(y), sep=""))
            ## if (!isCurrent(x)[["MSnSet"]])
            ##     x <- updateObject(x)
            assayData(x) <- combine(assayData(x), assayData(y))
            phenoData(x) <- combine(phenoData(x), phenoData(y))
            featureData(x) <- combine(featureData(x), featureData(y))
            experimentData(x) <- combine(experimentData(x),experimentData(y))
            protocolData(x) <- combine(protocolData(x), protocolData(y))
            x@processingData <- combine(processingData(x), processingData(y))
            x@processingData@processing <- paste("Combined [",
                                                 paste(dim(x), collapse = ","),
                                                 "] and [",
                                                 paste(dim(y), collapse = ","),
                                                 "] MSnSets ", date(), sep = "")
            x@qual <- data.frame() ## dropping qual slot
            ## annotation -- constant / not used
            if (validObject(x))
              return(x)
          })


setMethod("topN", signature(object = "matrix"), 
          function(object, groupBy, n=3, fun, ...) {
            if (missing(groupBy))
              stop("Specify how to group features to select top ", n, ".")
            if (missing(fun)) {
              fun <- sum
              if (ncol(object) > 1)
                message("Ranking features using their sum.")
            }
            rn <- rownames(object)
            idx <- by(object, groupBy, getTopIdx, n, fun, ...)
            object <- subsetBy(object, groupBy, idx)
            if (!is.null(rn)) {
              rownames(object) <- subsetBy(rn, groupBy, idx)
            } else {
              rownames(object) <- NULL
            }
            return(object)
          })

setMethod("topN", signature(object = "MSnSet"), 
          function(object, groupBy, n=3, fun, ...) {
            if (missing(groupBy))
              stop("Specify how to group features to select top ", n, ".")
            if (missing(fun)) {
              fun <- sum
              if (ncol(object) > 1)
                message("Ranking features using their sum.")
            }
            idx <- by(exprs(object), groupBy, getTopIdx, n, fun, ...)
            fn <- subsetBy(featureNames(object), groupBy, idx)
            .eset <- subsetBy(exprs(object), groupBy, idx)
            rownames(.eset) <- fn
            .proc <- processingData(object)
            .proc@processing <- c(.proc@processing,
                                  paste0("Selected top ", n,
                                         " features: ", date()))
            .fdata <- subsetBy(fData(object), groupBy, idx)
            message("Dropping spectrum-level 'qual' slot.")
            ans <- new("MSnSet",
                       experimentData = experimentData(object),
                       processingData = .proc,
                       exprs = .eset,
                       phenoData = phenoData(object),
                       featureData = new("AnnotatedDataFrame", data = .fdata),
                       annotation = object@annotation,
                       protocolData = protocolData(object))
            fn <- subsetBy(featureNames(object), groupBy, idx)
            if (validObject(ans))
              return(ans)
          })


setMethod("plotNA", signature(object = "MSnSet"), 
          function(object, pNA = .5) {
            if (pNA > 1)
              pNA <- 1
            if (pNA < 0)
              pNA <- 0
            X <- exprs(object)
            p <- plotNA.matrix(X, pNA)
            invisible(p)
          })


setMethod("plotNA", signature(object = "matrix"), 
          function(object, pNA = .5) {
            if (pNA > 1)
              pNA <- 1
            if (pNA < 0)
              pNA <- 0
            p <- plotNA.matrix(object, pNA)
            invisible(p)
          })

setMethod("filterNA", signature(object = "matrix"), 
          function(object, pNA = .5) {
            if (pNA > 1)
              pNA <- 1
            if (pNA < 0)
              pNA <- 0            
            k <- apply(object, 1,
                       function(x) sum(is.na(x))/length(x))
            accept <- k <= pNA
            if (sum(accept) == 1) {              
              ans <- matrix(object[accept, ], nrow=1)
              rownames(ans) <- rownames(object)[accept]
            } else {
              ans <- object[accept, ]
            }
            return(ans)
          })


setMethod("filterNA", signature(object = "MSnSet"), 
          function(object, pNA = .5) {
            if (pNA > 1)
              pNA <- 1
            if (pNA < 0)
              pNA <- 0
            k <- apply(exprs(object), 1,
                       function(x) sum(is.na(x))/length(x))
            accept <- k <= pNA
            object@processingData@processing <-
              c(processingData(object)@processing,
                paste0("Removed features with more that ",
                       pNA, "NAs: ", date()))
            ans <- object[accept, ]
            if (validObject(ans))
              return(ans)
          })

is.na.MSnSet <- function(x) is.na(exprs(x))
