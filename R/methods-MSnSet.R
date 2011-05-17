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
  msg <- validMsg(NULL, Biobase::assayDataValidMembers(assayData(object), c("exprs")))
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
              stop("Impurity matrix show be ",ncol(object)," by ",ncol(object))
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
setMethod("ratios","MSnSet",function(object,...) ratios.MSnSet(object,...))
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
