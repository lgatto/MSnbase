setMethod("all.equal", c("MSnExp", "MSnExp"),
          function(target, current, ...) {
              target@processingData@processing <-
                  current@processingData@processing <- NA_character_
              callNextMethod(target, current, ...)
          })

## The function below work for in memory and on disk MSnExp instances
## and for spectra lists.
equalMSnExps <- function(inmem, ondisk, ...) {
    if ((l <- length(inmem)) != length(ondisk))
        return("Objects have different lengths.")
    if (inherits(ondisk, "OnDiskMSnExp"))
        ondisk <- spectra(ondisk)
    if (inherits(inmem, "OnDiskMSnExp"))
        inmem <- spectra(inmem)
    ## otherwise, expected to be a list of spectra
    for (i in seq_len(l)) {
        if (!isTRUE(speq <- all.equal(ondisk[[i]], inmem[[i]], ...)))
            return(paste("Spectra at position", i, "are different."))
    }
    return(TRUE)
}

setMethod("all.equal", c("MSnExp", "OnDiskMSnExp"),
          function(target, current, ...) equalMSnExps(target, current, ...))

setMethod("all.equal", c("OnDiskMSnExp", "MSnExp"),
          function(target, current, ...) equalMSnExps(target, current, ...))

setMethod("all.equal", c("OnDiskMSnExp", "OnDiskMSnExp"),
          function(target, current, ...) {
              current@processingData@processing <-
                  target@processingData@processing <- NA_character_
              sp1 <- spectra(target)
              sp2 <- spectra(current)
              msg <- equalMSnExps(sp1, sp2, ...)
              ## call super to test other slots
              msg <- c(msg, callNextMethod(target, current))
              if (is.logical(msg) && all(msg)) TRUE
              else msg
          })
