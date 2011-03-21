## TODO getters and setters for MIAPE-MS data

##################################################################
## Methods for MIAPE class
setMethod("show","MIAPE",function(object) show.MIAPE(object))

setMethod("msInfo","MIAPE",
          function(object) msInfo.MIAPE(object))


setMethod("abstract","MIAPE",function(object) object@abstract)
setMethod("samples","MIAPE",function(object) object@samples)
setMethod("pubMedIds","MIAPE",function(object) object@pubMedIds)

setReplaceMethod("pubMedIds","MIAPE",function(object,value){
   object@pubMedIds = value
   object
})

setMethod("otherInfo","MIAPE",function(object) object@other)

setMethod("expinfo","MIAPE",
   function(object) {
      tmp <- c(object@name, object@lab, object@contact, object@title, object@url)
    names(tmp) <- c("name","lab","contact","title","url")
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

