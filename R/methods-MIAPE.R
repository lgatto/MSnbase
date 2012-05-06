## TODO getters and setters for MIAPE-MS data

##################################################################
## Methods for MIAPE class
setMethod("show","MIAPE",function(object) {
  show.MIAPE(object)
  invisible(NULL)
})

setMethod("msInfo","MIAPE",
          function(object) msInfo.MIAPE(object))


setMethod("abstract","MIAPE",function(object) object@abstract)
setMethod("samples","MIAPE",function(object) object@samples)
setMethod("pubMedIds","MIAPE",function(object) object@pubMedIds)

setReplaceMethod("pubMedIds","MIAPE",function(object,value){
   object@pubMedIds = value
   object
})

setMethod("email","MIAPE", function(object) object@email)
setMethod("title","MIAPE", function(object) object@title)
setMethod("ionSource","MIAPE", function(object) object@ionSource)
setMethod("analyser","MIAPE", function(object) object@analyser)
setMethod("detectorType","MIAPE", function(object) object@detectorType)

setReplaceMethod("pubMedIds","MIAPE",function(object, value){
   object@pubMedIds = value
   object
})


setMethod("otherInfo","MIAPE",function(object) object@other)

setMethod("expinfo","MIAPE",
   function(object) {
      tmp <- c(object@name, object@lab, object@contact,
               object@email, object@title, object@url)
    names(tmp) <- c("name","lab","contact",
                    "email", "title","url")
    return(tmp)
   }
)

setMethod("notes", signature(object="MIAPE"),
          function(object) object@other)

setReplaceMethod("notes", signature(object="MIAPE", value="list"),
                 function(object, value) {
                     object@other <- value
                     object
                 })

setReplaceMethod("notes", signature(object="MIAPE", value="character"),
                 function(object, value) {
                     object@other <- append(object@other, value)
                     object
                 })


## Adapted from Biobase::combine("MIAME", "MIAME") 
setMethod("combine",
          c("MIAPE", "MIAPE"),
          function(x, y, ...) {
            if (identical(x,y))
              return (x)
            for (sl in names(getSlots(class(x)))) {
              if (identical(slot(x, sl),slot(y, sl)))
                next
              slot(x, sl) <-
                switch(sl,
                       ## multiple elements possible
                       ## shared slots with MIAME
                       name = ,
                       lab = ,
                       contact = ,
                       title = ,
                       url = ,
                       pubMedIds = ,
                       samples = ,
                       hybridizations = ,
                       normControls = ,
                       preprocessing = ,
                       other = ,
                       ## MIAPE specific
                       dataStamp = ,
                       instrumentModel = ,
                       instrumentManufacturer = ,
                       instrumentCustomisations = ,
                       softwareName = ,
                       softwareVersion = ,
                       switchingCriteria = ,
                       isolationWidth = ,
                       parameterFile = ,
                       ionSource = ,
                       ionSourceDetails = ,
                       analyser = ,
                       analyserDetails = ,
                       collisionGas = ,
                       collisionPressure = ,
                       collisionEnergy = ,
                       detectorType = ,
                       detectorSensitivity = {
                         c(slot(x, sl), slot(y, sl))
                       },
                       ## just a single entry
                       abstract = {
                         paste(slot(x, sl), slot(y, sl), collapse = "\n")
                       },
                       .__classVersion__ = {
                         stop("'MIAPE' objects have different class version strings")
                       },
                       ## unknown
                       {
                         warning("\n  unknown or conflicting information in MIAPE field '",
                                 sl,"'; using information from first object 'x'")
                         slot(x, sl)
                       })
            }
            if (validObject(x))
              return(x)
          })
