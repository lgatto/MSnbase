#' @include hidden_aliases.R
NULL

setClass("BackendMzR", contains = "Backend")

#' @rdname hidden_aliases
setMethod("backendReadSpectra", "BackendMzR", function(object,
                                                       spectraData, ...) {
    file_f <- factor(spectraData$fileIdx, levels = unique(spectraData$fileIdx))
    fls <- object@files[as.integer(levels(file_f))]
    if (length(fls) == 1)
        .spectra_from_file_mzR(fls, spectraData)
    else
        unlist(mapply(fls, split(spectraData, file_f),
                      FUN = .spectra_from_file_mzR, SIMPLIFY = FALSE,
                      USE.NAMES = FALSE), recursive = FALSE)
})

#' @rdname hidden_aliases
setMethod(
    "backendWriteSpectra",
    "BackendMzR",
    function(object, spectra, spectraData) {
        if (isMSnbaseVerbose()) {
            message("Can not make changes to spectrum data persistent: ",
                    "'BackendMzR' is read-only.")
        }
    object
})

#' @rdname Backend
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
