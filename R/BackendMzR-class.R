#' @include hidden_aliases.R
NULL

setClass("BackendMzR", contains = "Backend")

#' @rdname hidden_aliases
setMethod("backendReadSpectra", "BackendMzR", function(object, file,
                                                       spectraData, ...,
                                                       BPPARAM=bpparam()) {
    spd_list <- split(spectraData, f = spectraData$fileIdx)
    if (length(spd_list) != length(file))
        stop("Number of files in 'spectraData' has to match length of 'file'")
    res <- bpmapply(FUN = .spectra_from_file_mzR, file, spd_list,
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
    names(res) <- NULL
    res <- unlist(res, recursive = FALSE)
    res[rownames(spectraData)]
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
