##################################################################
## Methods for MSnProcess class
setMethod("initialize","MSnProcess",          
          function(.Object,...) {
            .Object <- callNextMethod(.Object,...)
            .Object@MSnbaseVersion <- as.character(packageDescription("MSnbase",fields="Version"))
            if (validObject(.Object))
              return(.Object)
          })

setMethod("show","MSnProcess",function(object) show.MSnProcess(object))

setMethod("fileNames",
          signature(object="MSnProcess"),
          function(object) object@files)
