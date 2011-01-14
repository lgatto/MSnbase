##################################################################
## Methods for MSnSet class
## setMethod("initialize", "MSnSet",
##           function(.Object,...) {
##             .Object <- callNextMethod()
##             experimentData(.Object) <- new("MIAPE")
##             .Object
##           })

setMethod("show","MSnSet",
          function(object) show.MSnSet(object))

setMethod("normalise","MSnSet",
          function(object,method=c("sum","max"))
          normalise.MSnSet(object,match.arg(method))
          )
setMethod("dim","MSnSet",function(x) dim(exprs(x)))
setMethod("ratios","MSnSet",function(object) ratios.MSnSet(object))
setMethod("qual","MSnSet", function(object) object@qual)
## Not sure about these...
## setMethod("featureNames<-","MSnSet",
##           function(object,value="character") object@features <- value)
setReplaceMethod("featureNames",
                 signature(object="MSnSet",
                           value="character"),
                 function(object, value) {
                   object@features = value
                   if (validObject(object))
                     return(object)
                 })

setMethod("proteomicsData","MSnSet",function(object) object@proteomicsData)
setMethod("proteomicsData<-","MSnSet",
          function(object,value="MIAPE") object@proteomicsData <- value)

setMethod("fileNames",
          signature(object="MSnSet"),
          function(object) processingData(object)@files)

setReplaceMethod("proteomicsData",
                 signature(object="MSnSet",
                           value="MIAPE"),
                 function(object, value) {
                   object@proteomicsData = value
                   if (validObject(object))
                     return(object)
                 })


setMethod("processingData",
          signature(object="pSet"),
          function(object) object@process)
