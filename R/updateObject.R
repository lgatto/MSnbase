## Spectrum2: version 0.1.0 -> 0.2.0
## Change: new polarity slot inherited from Spectrum (previously in
##         Spectrum 1)
## Method: default updateObject works

setMethod("updateObject", signature(object = "MSnExp"),
          function(object, ..., verbose = FALSE) {
              if (verbose) message("updateObject(object = 'MSnExp')")              
              object <- asS4(object)
              ## if (isVersioned(object) && isCurrent(object)["MSnExp"])
              ##     return(callNextMethod())
              if (classVersion(object)["MSnExp"] == "0.3.0") {
                  e <- new.env()
                  for (sp in ls(assayData(object)))
                      assign(sp,
                             updateObject(get(sp, assayData(object))),
                             envir = e)
                  lockEnvironment(e, bindings = TRUE)
                  obj2 <- object
                  obj2@assayData <- e
                  obj2@.__classVersion__ <- classVersion(new("MSnExp"))
                  ## else if (classVersion(object)["MSnExp"] == "x.y.z") ...
              } else {
                  stop("Cannot update object of class '", class(object),
                       "', claiming to be MSnExp version '",
                       as(classVersion(object)["MSnExp"], "character"), "'")
              }
              if (validObject(obj2))
                  return(obj2)
          })
