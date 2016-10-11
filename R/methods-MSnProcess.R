##################################################################
## Methods for MSnProcess class
setMethod("initialize", "MSnProcess",
          function(.Object, ...) {
              .Object <- callNextMethod(.Object, ...)
              .Object@MSnbaseVersion <-
                  as.character(packageDescription("MSnbase",
                                                  fields = "Version"))
              if (validObject(.Object))
                  return(.Object)
          })

setMethod("show","MSnProcess",
          function(object) {
            show_MSnProcess(object)
            invisible(NULL)
          })

setMethod("fileNames",
          signature(object="MSnProcess"),
          function(object) object@files)

## setReplaceMethod("fileNames",
##           signature(object="MSnProcess", value="character"),
##           function(object, value) {
##             isExisiting <- file.exists(value)
##             if (!any(isExisiting))
##               stop("File(s) ", sQuote(value[!isExisiting]), " does not exist!")
##             object@files <- normalizePath(value)
##             return(object)
##           })

## Adapted from Biobase::combine("MIAME", "MIAME") 
setMethod("combine",
          c("MSnProcess", "MSnProcess"),
          function(x, y, ...) {
            if (identical(x,y))
              return (x)
            for (sl in names(getSlots(class(x)))) {
              if (identical(slot(x, sl),slot(y, sl)))
                next
              slot(x, sl) <-
                switch(sl,
                       files = ,
                       merged = ,
                       cleaned = ,
                       removedPeaks = ,
                       smoothed = ,
                       trimmed = ,
                       normalised = ,
                       MSnbaseVersion = {
                         c(slot(x, sl), slot(y, sl))
                       },
                       processing = {
                         paste("Combined MSnSets ", date(), sep = "")
                       },
                       .__classVersion__ = {
                         stop("'MSnProcess' objects have different class version strings")
                       },
                       {
                         warning("\n  unknown or conflicting information in MSnProcess field '",
                                 sl,"'; using information from object ", x)
                         slot(x, sl)
                       })
            }
            if (validObject(x))
              return(x)
          })
