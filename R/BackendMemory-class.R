#' @include hidden_aliases.R
NULL

#' BackendMemory class
#'
#' This class offers the *in-memory* reprensentation for the
#' [MSnExperiment-class] data. It mimics the classical [MSnExp-class] behaviour.
#'
#' @name BackendMemory-class
#' @docType class
#' @slot assay A `list` containing the [Spectrum-class] objects.
#' @family Backend classes
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @export
setClass("BackendMemory",
    contains="Backend",
    slots=c(
        assay="list" # or environment
    )
)

#' @describeIn BackendMemory-class Constructor
#'
#' This function is used to generated an *in-memory* backend. Just useful as
#' argument in [readMSnExperiment()].
#' @return A [BackendMemory-class].
#' @export
BackendMemory <- function() { new("BackendMemory") }

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendInitialize",
    signature="BackendMemory",
    definition=function(object, files, spectraData, ..., BPPARAM=bpparam()) {

    object@assay <- vector(mode="list", length=nrow(spectraData))
    names(object@assay) <-
        paste(.vdigest(files)[spectraData$fileIdx], spectraData$spIdx, sep="/")
    object
})

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendImportData",
    signature="BackendMemory",
    definition=function(object, files, spectraData, ..., BPPARAM=bpparam()) {

    spd <- split(spectraData, spectraData$fileIdx)

    split(object@assay, spectraData$fileIdx) <- bpmapply(
        .spectra_from_file_mzR, file=files, spectraData=spd,
        USE.NAMES=FALSE, SIMPLIFY=FALSE, BPPARAM=BPPARAM
    )
    object
})

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendReadSpectra",
    signature="BackendMemory",
    definition=function(object, file, spectraData, ...,
                        BPPARAM=bpparam()) {
    nms <- paste(.vdigest(file), spectraData$spIdx, sep="/")
    object@assay[nms]
})

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendWriteSpectra",
    signature="BackendMemory",
    definition=function(object, file, spectra, spectraData, ...,
                        BPPARAM=bpparam()) {
    nms <- paste(.vdigest(file), spectraData$spIdx, sep="/")
    object@assay[nms] <- spectra
    object
})
