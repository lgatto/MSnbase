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
    object@spectra <- lapply(
        object@spectra[rownames(spectraData)],
        function(s) {
            if (!is.null(s))
                s@fromFile <- match(s@fromFile, fidx)
            s
        }
    )
    callNextMethod()
})

#' @rdname hidden_aliases
setReplaceMethod(
    "backendSplitByFile",
    "BackendMemory",
    function(object, spectraData, ..., value) {
    rn <- split(rownames(spectraData), spectraData$fileIdx)
    fidx <- as.integer(sort.int(unique(spectraData$fileIdx)))
    for (i in seq(along=value)) {
        object@spectra[rn[[i]]] <-
            lapply(value[[i]]@spectra, function(s){ s@fromFile <- fidx[i]; s })
    }
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
    "backendReadSpectra",
    signature = "BackendMemory",
    definition = function(object, spectraData, ..., BPPARAM = bpparam()) {
    object@spectra[rownames(spectraData)]
})

#' @rdname hidden_aliases
setMethod(
    "backendWriteSpectra",
    signature = "BackendMemory",
    definition = function(object, spectra, spectraData, updateModCount, ...) {
        object@spectra[rownames(spectraData)] <- spectra
        if (updateModCount) {
            idx <- unique(vapply(spectra, fromFile, integer(1L)))
            object@modCount[idx] <- object@modCount[idx] + 1L
        }
        validObject(object)
        object
})

#' @rdname hidden_aliases
setMethod("backendUpdateMetadata", "BackendMemory", function(object,
                                                             spectraData) {
    object@spectra <- object@spectra[rownames(spectraData)]
    object@spectra <- mapply(object@spectra,
                             split(spectraData, seq_len(nrow(spectraData))),
                             FUN = .spectrum_set_header)
    validObject(object)
    object
})
