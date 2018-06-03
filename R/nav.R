## TODO: - make it work with an Spectrum object, so that one can chain
## the calls.

##' Navigate an \code{MSnExp} object by moving to the next or previous
##' spectrum.
##'
##' @title Navigate an \code{MSnExp} object
##' @param i The name or index of the current spectrum
##' @param object The \code{MSnExp} object
##' @param msLevel The MS level of the next or previous spectrum. If
##'     missing (default), then the level of the current spectrum is
##'     used.
##' @param nav One of \code{"nextMS"} or \code{"prevMS"}, to obtain
##'     the next or previous spectrum of level \code{msLevel}.
##' @param ... Additional parameters. Currently ignored.
##' @return An object of class \code{Spectrum1} or \code{Spectrum2},
##'     depending on the value of \code{msLevel} or \code{NULL}, of no
##'     spectrum is found.
##' @author Laurent Gatto
##' @examples
##' f <- msdata::proteomics(full.names = TRUE, pattern = "MS3")
##' x <- readMSData(f, centroided. = c(FALSE, TRUE, FALSE), mode = "onDisk")
##' (sp <- which(msLevel(x) == 3)[2]) ## 2nd MS3 spectrum
##' x[[sp]] ## curent MS3
##' MSnbase:::nextMS(sp, x) ## next MS3
##' MSnbase:::prevMS(sp, x) ## prev MS3
##' MSnbase:::prevMS(sp, x, 2L) ## prev MS2
##' MSnbase:::prevMS(sp, x, 1L) ## prev MS1
navMS <- function(i,
                  object,
                  msLevel,
                  nav = c("nextMS", "prevMS"),
                  ...) {
    nav <- match.arg(nav)
    if (is.character(i))
        i <- which(featureNames(object) == i)
    stopifnot(length(i) == 1 && is.numeric(i))
    if (missing(msLevel))
        msLevel <- fData(object)[i, "msLevel"]
    ks <- which(fData(object)[, "msLevel"] == msLevel)
    k <- switch(nav,
                nextMS = head(ks[ks > i], n = 1),
                prevMS = tail(ks[ks < i], n = 1))
    if (length(k))
        return(object[[k]])
    return(NULL)
}


##' @rdname navMS
nextMS <- function(...) navMS(..., nav = "nextMS")

##' @rdname navMS
prevMS <- function(...) navMS(..., nav = "prevMS")
