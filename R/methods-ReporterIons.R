##################################################################
## Methods for ReporterIons class
setMethod("show","ReporterIons",
          function(object) {
            show_ReporterIons(object)
            invisible(NULL)
          })

setMethod("[","ReporterIons",
          function(x,i,j="missing",drop="missing") "[.ReporterIons"(x,i))
setMethod("length","ReporterIons",function(x) length(x@mz))
setMethod("mz","ReporterIons", function(object) object@mz)

setMethod("width","ReporterIons", function(x) x@width)

setMethod("reporterColours","ReporterIons", function(object) object@col)
setMethod("reporterColors" ,"ReporterIons", function(object) object@col)

setMethod("reporterNames","ReporterIons", function(object) object@reporterNames)

setReplaceMethod("reporterNames",
                 signature(object="ReporterIons",
                           value="character"),
                 function(object, value) {
                   if (length(value)!=length(object))
                     stop(paste("Please provide names for",
                                length(object),
                                "reporters",sep=" "))
                   object@reporterNames = value
                   if (validObject(object))
                     return(object)
                 })

setMethod("initialize", "ReporterIons",
          function(.Object,...) {
            .Object <- callNextMethod()
            if (length(.Object@mz)!=length(.Object@reporterNames)) {
              warning("Setting reporter name.")
              .Object@reporterNames <- paste(.Object@name,.Object@mz,sep=".")
            }
            if (validObject(.Object))
              return(.Object)
          })


setMethod("description", "ReporterIons",
          function(object) object@description)

setMethod("names", "ReporterIons",
          function(x) x@name)
