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

.valid.BackendMemory.match.file.spectra <- function(files, spectra) {
    n <- c(length(files), length(spectra))

    if (all(n)) {
        fnms <- names(files)
        snms <- names(spectra)

        if (isFALSE(all(.BackendMemory.fileIndexFromName(snms) %in% fnms))) {
            return("Mismatch between names for 'files' and 'spectra' found.")
        }
    }
    NULL
}

setValidity("BackendMemory", function(object) {
    lapply(object@spectra, validObject)
    msg <- c(
        .valid.BackendMemory.spectra.names(object@spectra)
        ## .valid.BackendMemory.match.file.spectra(object@files, object@spectra)
    )

    if (is.null(msg)) { TRUE } else { msg }
})

#' @rdname Backend
BackendMemory <- function() { new("BackendMemory") }

#' @rdname hidden_aliases
setMethod("backendSubset", "BackendMemory", function(object, spectraData) {
    files_orig <- object@files
    object@files <- object@files[unique(spectraData$fileIdx)]
    ## Update also `@fromFile` in the spectra.
    object@spectra <- lapply(object@spectra[rownames(spectraData)],
                             function(z) {
                                 z@fromFile <- match(files_orig[z@fromFile],
                                                     object@files)
                                 z
                             })
    validObject(object)
    object
})

## #' @rdname hidden_aliases
## setMethod(
##     "[",
##     signature("BackendMemory", i="ANY", j="missing"),
##     function(x, i, j, ..., drop=FALSE) {
##     x@spectra <- x@spectra[i]
##     f <- .BackendMemory.fileIndexFromName(names(x@spectra))
##     x@files <- x@files[f]
##     validObject(x)
##     x
## })

## #' @rdname hidden_aliases
## setMethod("filterFile", "BackendMemory", function(object, file, ...) {
##     ## we don't use `callNextMethod` here, because it will double the
##     ## `validObject` call and the first `validObject` fails with an error
##     ## because the names of `file` and `spectra` don't match anymore (`file` is
##     ## already filtered but `spectra` is not)
##     if (is.character(file)) {
##         file <- base::match(file, object@files)
##     }
##     object@files <- object@files[file]
##     keep <- .BackendMemory.fileIndexFromName(names(object@spectra)) %in%
##         names(object@files)
##     object@spectra <- object@spectra[keep]
##     validObject(object)
##     object
## })

#' @rdname hidden_aliases
setMethod(
    "backendInitialize",
    signature="BackendMemory",
    definition=function(object, files, spectraData, ..., BPPARAM=bpparam()) {

    object@spectra <- vector(mode="list", length=nrow(spectraData))
    names(object@spectra) <- rownames(spectraData)
    object <- callNextMethod()
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "backendImportData",
    signature="BackendMemory",
    definition=function(object, spectraData, ..., BPPARAM=bpparam()) {

    spd <- split(spectraData, spectraData$fileIdx)

    split(object@spectra, spectraData$fileIdx) <- bpmapply(
        .spectra_from_file_mzR, file=object@files, spectraData=spd,
        USE.NAMES=FALSE, SIMPLIFY=FALSE, BPPARAM=BPPARAM
    )
    validObject(object)
    object
})

#' @rdname hidden_aliases
setMethod(
    "backendReadSpectra",
    signature="BackendMemory",
    definition=function(object, spectraData, ...,
                        BPPARAM=bpparam()) {
        object@spectra[rownames(spectraData)]
    })

#' @rdname hidden_aliases
setMethod(
    "backendWriteSpectra",
    signature="BackendMemory",
    definition=function(object, spectra, spectraData, ...,
                        BPPARAM=bpparam()) {
        object@spectra[rownames(spectraData)] <- spectra
        validObject(object)
        object
})

.BackendMemory.fileIndexFromName <- function(x) {
    gsub("\\.S[0-9]+$", "", x)
}
