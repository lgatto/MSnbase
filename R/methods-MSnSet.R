##################################################################
## Methods for MSnSet class
setMethod("initialize", "MSnSet",
          function(.Object,experimentData,...) {
            .Object <- callNextMethod()
            if (missing(experimentData)) 
              experimentData(.Object) <- new("MIAPE")
            else
              experimentData(.Object) <- experimentData
            .Object
          })

                   
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
setMethod("featureNames<-","MSnSet",
          function(object,value="character") object@features <- value)
setReplaceMethod("featureNames",
                 signature(object="MSnSet",
                           value="character"),
                 function(object, value) {
                   object@features = value
                   if (validObject(object))
                     return(object)
                 })


