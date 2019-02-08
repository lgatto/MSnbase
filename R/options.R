##' MSnbase defined a few options globally using the standard R
##' options mechanism. The current values of these options can be
##' queried with `MSnbaseOptions`. The options are:
##' + `verbose`: defines a session-wide verbosity flag, that
##'   is used if the `verbose` argument in individual functions is
##'   not set.
##' + `PARALLEL_THRESH`: defines the minimum number of spectra per file
##'   necessary before using parallel processing.
##' + `fastLoad`: `logical(1)`. If `TRUE` performs faster data loading for all
##'   methods of [OnDiskMSnExp] that load data from the original files (such as
##'   [spectrapply()]). Users experiencing data I/O errors (observed mostly
##'   on macOS systems) should set this option to `FALSE`.
##' + `HDF5_COMP_LEVEL`: defines the compression level of the HDF5 files.
##'   Supports integer values between 0 (no compression) and 9 (highest
##'   compression).
##'
##' `isMSnbaseVerbose` is one wrapper for the verbosity flag,
##' also available through `options("MSnbase")$verbose`.
##'
##' There are also setters to set options individually. When run
##' without argument, the verbosity setter inverts the current value
##' of the option.
##'
##' @title MSnbase options
##'
##' @param opt The value of the new option
##'
##' @return A `list` of MSnbase options and the single option
##'     values for the individual accessors.
##' @md
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

##' @rdname MSnbaseOptions
setMSnbaseFastLoad <- function(opt = TRUE) {
    oldopts <- opts <- MSnbaseOptions()
    opts$fastLoad <- opt
    options("MSnbase" = opts)
    invisible(oldopts$fastLoad)
}

##' @rdname MSnbaseOptions
isMSnbaseFastLoad <- function() {
    fast_load <- MSnbaseOptions()$fastLoad
    ## For some odd reasons we get also NULL back - parallel processing?
    if (!length(fast_load))
        fast_load <- FALSE
    fast_load
}

.hdf5_compression_level <- function() {
    MSnbaseOptions()$HDF5_COMP_LEVEL
}

##' @rdname MSnbaseOptions
setHdf5CompressionLevel <- function(opt) {
    opt <- as.integer(opt)
    if (!(opt %in% 0:9))
        stop("Compression level should be an integer between 0 and 9")
    opts <- MSnbaseOptions()
    opts$HDF5_COMP_LEVEL <- opt
    options("MSnbase" = opts)
    invisible(opts$HDF5_COMP_LEVEL)
}
