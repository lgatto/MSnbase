##' MSnbase defined a few options globally using the standard R
##' options mechanism. The current values of these options can be
##' queried with \code{MSnbaseOptions}. The options are
##' \code{verbose}, that defines a session-wide verbosity flag, that
##' is used if the \code{verbose} argument in individual functions is
##' not set. The \code{PARALLEL_THRESH} defines the minimum number of
##' spectra per file necessary before using parallel processing.
##'
##' \code{isMSnbaseVerbose} is one wrapper for the verbosity flag,
##' also available through \code{options("MSnbase")$verbose}.
##'
##' There are also setters to set options individually. When run
##' without argument, the verbosity setter inverts the current value
##' of the option.
##' 
##' @title MSnbase options
##' @param opt The value of the new option
##' @return A \code{list} of MSnbase options and the single option
##'     values for the individual accessors.
MSnbaseOptions <- function()
    options("MSnbase")[[1]]

##' @rdname MSnbaseOptions
isMSnbaseVerbose <- function()
    MSnbaseOptions()$verbose

##' @rdname MSnbaseOptions
setMSnbaseVerbose <- function(opt) {
    if (missing(opt))
        opt <- !isMSnbaseVerbose()
    oldopts <- opts <- MSnbaseOptions()
    opts$verbose <- opt
    options("MSnbase" = opts)
    invisible(oldopts$verbose)
}

##' @rdname MSnbaseOptions
setMSnbaseParallelThresh <- function(opt = 1000) {
    oldopts <- opts <- MSnbaseOptions()
    opts$PARALLEL_THRESH <- opt
    options("MSnbase" = opts)
    invisible(oldopts$PARALLEL_THRESH)
}
