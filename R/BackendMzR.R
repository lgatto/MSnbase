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
