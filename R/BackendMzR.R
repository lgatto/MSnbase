#' @include hidden_aliases.R
NULL

#' @title Mass spectrometry data files-based backend
#'
#' @description
#'
#' The `BackendMzR` uses the original MS data files (such as *mzML*, *mzXML* or
#' *CDF* files) as backend and reads the data on demand from these files. This
#' ensures a low memory footprint and enables thus the analysis also of very
#' large experiments - at the cost of a slightly lower performance.
#'
#' @name BackendMzR-class
#'
#' @author Johannes Rainer
#'
#' @family backends
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' ## TODO: define to create the documentation. Add this to the general
#' ## backend documentation or have it in its own one.
setClass("BackendMzR", contains = "Backend")

setMethod("backendReadSpectra", "BackendMzR", function(object, files,
                                                       spectraData, ...) {
})

#' @rdname BackendMzR-class
#'
#' @md
#'
#' @noRd
BackendMzR <- function() {
    new("BackendMzR")
}

#' @description
#'
#' Internal function to read spectrum data from MS raw files (using mzR) and
#' return a list of `Spectrum` objects.
#'
#' @param file `character(1)` with the (absolute path) to the MS raw data file.
#'
#' @param spectraData `data.frame` or `DataFrame` with the spectrum data for
#'     spectra **of a single MS file**.
#'
#' @return `list` of `Spectrum` objects.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.spectra_from_file_mzR <- function(file, spectraData) {
    if (missing(file) || missing(spectraData))
        stop("Both 'file' and 'spectraData' are required")
    fl <- .openMSfile(file)
    on.exit(mzR::close(fl))
    hdr <- mzR::header(fl, max(spectraData$spIdx))
    mzi <- mzR::peaks(fl, spectraData$spIdx)
    if (is.matrix(mzi))
        mzi <- list(mzi)
    .spectra_from_data(mzi, spectraData)
}
