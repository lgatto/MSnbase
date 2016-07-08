.slotNames <- function(x, toremove = NULL) {
    toremove <- c(toremove, ".__classVersion__")
    setdiff(slotNames(x), toremove)
}

setMethod("updateObject", signature(object = "Spectrum"),
          function(object, ..., verbose=FALSE) {
              if (verbose) message("updateObject(object = 'Spectrum')")
              object <- asS4(object)
              if (isVersioned(object) && isCurrent(object)["Spectrum"])
                callNextMethod()
              else {
                  ## version 0.3.0. -> 0.4.0 needs a new slot "smoothed"
                  to <- new(class(object))
                  for (sl in .slotNames(object, "smoothed"))
                      slot(to, sl) <- slot(object, sl)
                  if (validObject(to))
                      to
              }
          })


setMethod("updateObject", signature(object = "MSnExp"),
          function(object, ..., verbose = FALSE) {
              if (verbose) message("updateObject(object = 'MSnExp')")
              object <- asS4(object)
              if (isVersioned(object) && isCurrent(object)["MSnExp"])
                  callNextMethod()
              else {
                  to <- new(class(object))
                  for (sl in .slotNames(object))
                      slot(to, sl) <- slot(object, sl)
                  e <- list2env(lapply(spectra(object), updateObject))
                  lockEnvironment(e, bindings = TRUE)
                  to@assayData <- e
                  to <- logging(to,
                                paste("Updated from version",
                                      as(classVersion(object)["MSnExp"], "character"),
                                      "to",
                                      as(classVersion(to)["MSnExp"], "character")))
                  if (validObject(to))
                      to
              }
          })
