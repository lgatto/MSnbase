setMethod("all.equal", c("MSnExp", "MSnExp"),
          function(target, current, ...) {
              target@processingData@processing <- 
                  current@processingData@processing <- NA_character_
              callNextMethod(target, current)
          })

setMethod("all.equal", c("MSnExp", "OnDiskMSnExp"),
          function(target, current, ...) equalMSnExps(target, current))
setMethod("all.equal", c("OnDiskMSnExp", "MSnExp"),
          function(target, current, ...) equalMSnExps(current, target))

equalMSnExps <- function(inmem, ondisk, ...) {
    if ((l <- length(inmem)) != length(ondisk)) 
        return("Objects have different lengths.")
    if (!isTRUE(fneq <- all.equal(featureNames(inmem), featureNames(ondisk)))) 
        return("Object have different feature names.")
    ## access data on disk only once
    spondisk <- spectra(ondisk)
    for (i in seq_len(l)) {
        if (!isTRUE(speq <- all.equal(spondisk[[i]], inmem[[i]]))) 
            return("Spectra at position", i, "are different.")
    }
    return(TRUE)
}


setMethod("all.equal", c("OnDiskMSnExp", "OnDiskMSnExp"),
          function(target, current, ...) {
              current@processingData@processing <-
                  target@processingData@processing <- NA_character_
              sp1 <- spectra(target)
              sp2 <- spectra(current)
              msg <- all.equal(sp1, sp2)
              msg <- c(msg, callNextMethod(target, current))
              if (is.null(msg)) TRUE
              else msg 
          })




