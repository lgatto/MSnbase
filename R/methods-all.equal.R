setMethod("all.equal", c("MSnExp", "MSnExp"),
          function(target, current, ...) {
              target@processingData <- new("MSnProcess")
              current@processingData <- new("MSnProcess")
              callNextMethod(target, current)
          })

setMethod("all.equal", c("MSnExp", "OnDiskMSnExp"),
          function(target, current, ...) equalMSnExps(target, current))
setMethod("all.equal", c("OnDiskMSnExp", "MSnExp"),
          function(target, current, ...) equalMSnExps(current, target))

equalMSnExps <- function(inmem, ondisk, ...) {
    msg <- NULL
    if ((l <- length(inmem)) != length(ondisk)) {
        msg <- c(msg, "Objects of different length")
    } else {
        ## access data on disk only once
        spondisk <- spectra(ondisk)
        for (i in seq_len(l)) {
            if (!isTRUE(speq <- all.equal(spondisk[[i]], inmem[[i]]))) {
                msg <- c(msg, speq)
                break
            }
        }
    }
    if (is.null(msg)) TRUE
    else msg
}
