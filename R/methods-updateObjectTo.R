.slotNames0 <- function(x, toremove = NULL) {
    toremove <- c(toremove, ".__classVersion__")
    setdiff(slotNames(x), toremove)
}

## Spectrum2:
## Version 0.1.0 -> 0.2.0
## Change: new polarity slot inherited from Spectrum (previously in
##         Spectrum 1)
## Version 0.3.0. -> 0.4.0
## Change: Spectrum2 get new slot "smoothed" slot

setMethod("updateObject", signature(object = "Spectrum"),
          function(object, ..., verbose = isMSnbaseVerbose()) {
              if (verbose) message("updateObject(object = 'Spectrum')")
              object <- asS4(object)
              if (isVersioned(object) && isCurrent(object)["Spectrum"])
                callNextMethod()
              else {
                  to <- new(class(object))
                  if (object@.__classVersion__["Spectrum2"] == "0.1.0") {
                      torm <- c("smoothed", "polarity")
                  } else if (object@.__classVersion__["Spectrum2"] == "0.2.0") {
                      torm <- "smoothed"
                  }
                  for (sl in .slotNames0(object, torm))
                      slot(to, sl) <- slot(object, sl)
                  if (validObject(to))
                      to
              }
          })


setMethod("updateObject", signature(object = "MSnExp"),
          function(object, ..., verbose = isMSnbaseVerbose()) {
              if (verbose) message("updateObject(object = 'MSnExp')")
              object <- asS4(object)
              if (isVersioned(object) && isCurrent(object)["MSnExp"])
                  callNextMethod()
              else {
                  to <- new(class(object))
                  for (sl in .slotNames0(object))
                      slot(to, sl) <- slot(object, sl)
                  e <- list2env(lapply(spectra(object), updateObject))
                  lockEnvironment(e, bindings = TRUE)
                  to@assayData <- e
                  fromv <- as(classVersion(object)["MSnExp"], "character")
                  tov <- as(classVersion(to)["MSnExp"], "character")
                  to <- logging(to,
                                paste("Updated from version", fromv,
                                      "to", tov))
                  if (validObject(to))
                      to
              }
          })
