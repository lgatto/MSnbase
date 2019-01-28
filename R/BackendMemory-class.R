#' @include hidden_aliases.R
NULL

setClass("BackendMemory",
    contains="Backend",
    slots=c(
        spectra="list" # or environment
    )
)

.valid.BackendMemory.spectra.names <- function(x) {
    n <- length(x)
    nms <- names(x)

    if (n) {
        if (any(is.null(nms)))
            return("Spectra names should not be NULL.")
        if (anyNA(nms))
            return("Spectra names should not contain NA.")
        if (!all(nchar(nms)))
            return("Spectra names should not be missing.")
        if (anyDuplicated(nms))
            return("Duplicated spectra names found.")
    }
    NULL
}

setValidity("BackendMemory", function(object) {
    lapply(object@spectra, validObject)
    msg <- .valid.BackendMemory.spectra.names(object@spectra)

    if (is.null(msg)) { TRUE } else { msg }
})

#' This function is used to generated an `BackendMemory` object
#' (*in-memory* backend). It doesn't support any arguments.
#'
#' @return A [BackendMemory-class].
#' @rdname Backend
BackendMemory <- function() { new("BackendMemory") }

#' @rdname hidden_aliases
setMethod(
    "backendInitialize",
    signature="BackendMemory",
    definition=function(object, files, spectraData, ..., BPPARAM=bpparam()) {

    object@spectra <- vector(mode="list", length=nrow(spectraData))
    names(object@spectra) <-
        paste(.vdigest(files)[spectraData$fileIdx], spectraData$spIdx, sep="/")
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "backendImportData",
    signature="BackendMemory",
    definition=function(object, files, spectraData, ..., BPPARAM=bpparam()) {

    spd <- split(spectraData, spectraData$fileIdx)

    split(object@spectra, spectraData$fileIdx) <- bpmapply(
        .spectra_from_file_mzR, file=files, spectraData=spd,
        USE.NAMES=FALSE, SIMPLIFY=FALSE, BPPARAM=BPPARAM
    )
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "backendReadSpectra",
    signature="BackendMemory",
    definition=function(object, file, spectraData, ...,
                        BPPARAM=bpparam()) {
    nms <- paste(.vdigest(file), spectraData$spIdx, sep="/")
    object@spectra[nms]
})

#' @rdname hidden_aliases
setMethod(
    "backendWriteSpectra",
    signature="BackendMemory",
    definition=function(object, file, spectra, spectraData, ...,
                        BPPARAM=bpparam()) {
    nms <- paste(.vdigest(file), spectraData$spIdx, sep="/")
    object@spectra[nms] <- spectra
    validObject(object)
    object
})
