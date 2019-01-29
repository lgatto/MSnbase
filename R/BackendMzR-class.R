#' @include hidden_aliases.R
NULL

setClass("BackendMzR", contains = "Backend")

#' @rdname hidden_aliases
setMethod("backendSpectrapply", "BackendMzR", function(object, spectraData,
                                                       FUN = NULL, ...,
                                                       BPPARAM = bpparam()) {
    file_f <- factor(spectraData$fileIdx, levels = unique(spectraData$fileIdx))
    fls <- object@files[as.integer(levels(file_f))]
    if (any(is.na(fls)))
        stop("file index of 'spectraData' out of bounds")
    pqueue <- object@processingQueue
    if (!is.null(FUN))
        pqueue <- c(pqueue, list(ProcessingStep(FUN, ARGS = list(...))))
    unlist(bpmapply(split(spectraData, f = file_f), fls,
                    FUN = function(sp, fl, queue) {
                        .apply_processing_queue(.spectra_from_file_mzR(fl, sp),
                                                queue)
                    }, MoreArgs = list(queue = pqueue), BPPARAM = BPPARAM,
                    SIMPLIFY = FALSE, USE.NAMES = FALSE), recursive = FALSE)
})

#' @rdname hidden_aliases
setMethod("backendReadSpectra", "BackendMzR", function(object,
                                                       spectraData, ...,
                                                       BPPARAM=bpparam()) {
    backendSpectrapply(object, spectraData = spectraData, ...,
                       BPPARAM = BPPARAM)
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
