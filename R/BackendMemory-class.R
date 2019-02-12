#' @include hidden_aliases.R
NULL

setClass("BackendMemory",
    contains = "Backend",
    slots = c(
        spectra = "list" # or environment
    )
)

.valid.BackendMemory.spectra.names <- function(x) {
    n <- length(x)

    if (n) {
        nms <- names(x)

        if (any(is.null(nms)))
            return("Spectra names should not be NULL.")
        if (anyNA(nms))
            return("Spectra names should not contain NA.")
        if (!all(nchar(nms)))
            return("Spectra names should not be missing.")
        if (anyDuplicated(nms))
            return("Duplicated spectra names found.")
        if (isFALSE(all(grepl("^F[0-9]+\\.S[0-9]+$", nms))))
            return("Names of 'spectra' don't follow F[0-9]+.S[0-9]+ format.")
    }
    NULL
}

setValidity("BackendMemory", function(object) {
    lapply(object@spectra, validObject)
    msg <- .valid.BackendMemory.spectra.names(object@spectra)
    if (is.null(msg)) { TRUE } else { msg }
})

#' @rdname Backend
BackendMemory <- function() { new("BackendMemory") }

#' @rdname hidden_aliases
setMethod("backendSubset", "BackendMemory", function(object, spectraData) {
    fidx <- unique(spectraData$fileIdx)
    ## Update also `@fromFile` in the spectra.
    object@spectra <- lapply(object@spectra[rownames(spectraData)], function(z) {
        z@fromFile <- match(z@fromFile, fidx)
        z
    })
    callNextMethod()
})

#' @rdname hidden_aliases
setMethod(
    "backendInitialize",
    signature = "BackendMemory",
    definition = function(object, files, spectraData, ..., BPPARAM = bpparam()) {

    object@spectra <- vector(mode = "list", length = nrow(spectraData))
    names(object@spectra) <- rownames(spectraData)
    object <- callNextMethod()
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "backendImportData",
    signature = "BackendMemory",
    definition = function(object, spectraData, ..., BPPARAM = bpparam()) {

    spd <- split(spectraData, spectraData$fileIdx)

    split(object@spectra, spectraData$fileIdx) <- bpmapply(
        .spectra_from_file_mzR, file = object@files, spectraData = spd,
        USE.NAMES = FALSE, SIMPLIFY = FALSE, BPPARAM = BPPARAM
    )
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "backendReadSpectra",
    signature = "BackendMemory",
    definition = function(object, spectraData, ..., BPPARAM = bpparam()) {
    object@spectra[rownames(spectraData)]
})

#' @rdname hidden_aliases
setMethod(
    "backendWriteSpectra",
    signature = "BackendMemory",
    definition = function(object, spectra, spectraData, ...,
                          BPPARAM = bpparam()) {
    object@spectra[rownames(spectraData)] <- spectra
    validObject(object)
    object
})
