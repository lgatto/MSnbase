#' @include hidden_aliases.R
NULL

#' @title Mass spectrometry data managing backends
#'
#' @aliases Backend-class BackendMzR-class BackendMemory-class BackendHdf5-class
#'
#' @description
#'
#' [MSnExperiment-class] objects support the use of different backends to
#' manage and access mass spectrometry data. Backends can be generally
#' classified into *in-memory* and *on-disk* backends. In-memory backends keep
#' all the (spectra) data in memory ensuring fast data access. On-disk backends
#' do not keep any data in memory but fetch the requested spectrum data only on
#' demand, applying eventual data manipulations on-the-fly. Due to their minimal
#' memory demand, on-disk backends support also loading and analyzing very large
#' MS experiments.
#'
#' Available backends in `MSnbase` are listed below.
#'
#' @section BackendMemory:
#'
#' The `BackendMemory` uses a `list` as backend and stores all MS data in the
#' memory. This ensures a high performance but needs a lot of memory for larger
#' experiments.  It mimics the classical [MSnExp-class] behaviour.
#'
#' New backends can be created with the `BackendMemory()` function.
#
#' @section BackendMzR:
#'
#' The `BackendMzR` uses the original MS data files (such as *mzML*, *mzXML* or
#' *CDF* files) as backend and reads the data on demand from these files. This
#' ensures a low memory footprint and enables thus the analysis also of very
#' large experiments - at the cost of a slightly lower performance. New
#' backends can be created with the `BackendMzR` function.
#'
#' New backends can be created with the `BackendMzR()` function.
#'
#' @section BackendHdf5:
#'
#' The `BackendHdf5` is, similar to the `BackendMzR`, a *on-disk* backend that
#' does only keep the minimum required data in memory (i.e. spectrum metadata).
#' The m/z and intensity values of all spectra are stored in HDF5 files (one
#' for each input file). This backend combines the advantages of the
#' `BackendMzR` (low memory footprint) with faster data access and the support
#' to apply data manipulations persistently to the data. Also, reading data
#' from HDF5 files is considerably faster than reading data from MS raw files
#' (mzML, mzXML or CDF). By default, HDF5 files are stored in the current
#' working directory, but it is also possible to specify a directory with
#' the `path` parameter of the [readMSnExperiment()] function (passed as an
#' optional parameter).
#'
#' New backends can be created with the `BackendHdf5()` function.
#'
#' @name Backend
#'
#' @author Sebastian Gibb, Johannes Rainer
#'
#' @md
NULL

#' Base class for all other [MSnExperiment-class] data
#' backend classes (e.g., [BackendMemory-class], [BackendHdf5-class]).
#'
#' It is a *VIRTUAL* class and just defines a common interface for possible
#' backends.
#'
#' @section Implementation notes:
#' An implementation *MUST* provide methods for the following generics:
#'
#' - [backendReadSpectra()] to read spectra from the backend.
#' - [backendWriteSpectra()] to write spectra to the backend.
#' - [backendSubset()] to subset the backend.
#'
#' It may also provide methods for:
#'
#' - [backendInitialize()] to setup the backend (create files, tables, ...).
#'
#' @name Backend
#'
#' @author Sebastian Gibb
#'
#' @noRd
setClass("Backend",
    slots = c(
        files = "character",    # src files (i.e. mzML files)
        # is incremented every time `backendWriteSpectra` is called and is used
        # to test for changes of the hdf5 files and superficial copying
        # see issue https://github.com/lgatto/MSnbase/issues/429
        modCount = "integer"
    ),
    prototype = prototype(files = character(), modCount = integer()),
    contains = "VIRTUAL"
)

.valid.Backend.files <- function(x) {
    n <- length(x)

    if (n) {
        if (anyNA(x))
            return("Files should not contain NA.")
        if (!all(nchar(x)))
            return("Files should not be missing.")
        if (anyDuplicated(x))
            return("Duplicated file names found.")
    }
    NULL
}

.valid.Backend.modCount <- function(x, y) {
    if (length(x) != length(y))
        "Different number of source files and modification counters."
    else
        NULL
}

setValidity("Backend", function(object) {
    msg <- c(.valid.Backend.files(object@files),
             .valid.Backend.modCount(object@files, object@modCount))
    if (is.null(msg)) { TRUE } else { msg }
})

#' @rdname hidden_aliases
#' @param object Object to display.
setMethod("show", signature = "Backend", definition = function(object) {
    cat("Backend:", class(object)[1L], "\n")
    fls <- basename(object@files)
    if (length(fls) > 3)
        fls <- c(fls[1:3], paste0("(", length(fls) - 3,
                                  " more. Use `fileNames` to list all.)"))
    cat("Source files:\n", paste(" ", fls, collapse = "\n"), "\n", sep = "")
})

#' @rdname hidden_aliases
setMethod("fileNames", "Backend", function(object, ...) object@files)

#' Initialize a backend
#'
#' This generic is used to setup a backend.
#'
#' It should be only reimplemented if the backend needs specific requirements,
#' e.g. for *HDF5* the creation of .h5 files; for *SQL* the creation of a
#' database or tables.
#'
#' @param object An object inheriting from [Backend-class],
#' i.e. [BackendHdf5-class]
#' @param files The path to the source (generally .mzML) files.
#' @param spectraData A [S4Vectors::DataFrame-class]
#' @param ... Other arguments passed to the methods.
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setGeneric(
    "backendInitialize",
    def = function(object, files, spectraData, ...)
        standardGeneric("backendInitialize"),
    valueClass = "Backend"
)
setMethod(
    "backendInitialize",
    signature = "Backend",
    definition = function(object, files, spectraData, ...) {
    object@files <- files
    object@modCount <- integer(length(files))
    validObject(object)
    object
})

#' Read spectra data from a backend
#'
#' This generic is used to read spectra data from a backend.
#'
#' It *MUST* be reimplemented by all backends!
#'
#' @inheritParams backendInitialize
#' @param file The path to the source (generally .mzML) file.
#' @param processingQueue `list` of `ProcessingStep` objects defining the
#'     processing steps to be applied to the spectra before returning them.
#' @return A list of [Spectrum-class] objects.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setGeneric(
    "backendReadSpectra",
    def = function(object, spectraData, processingQueue = list(), ...)
        standardGeneric("backendReadSpectra"),
    valueClass = "list"
)

#' Write spectra data to a backend
#'
#' This generic is used to write spectra data to a backend.
#'
#' It *MUST* be reimplemented by all backends!
#'
#' @inheritParams backendReadSpectra
#' @param spectra A `list` of [Spectrum-class] objects that should be written to
#' the backend.
#' @param updateModCount `logical`, should the `@modCount` be incremented? Only
#' set to `FALSE` if you know what you are doing (e.g. in `setBackend`).
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setGeneric(
    "backendWriteSpectra",
    def = function(object, spectra, spectraData, updateModCount = TRUE, ...)
        standardGeneric("backendWriteSpectra"),
    valueClass = "Backend"
)

#' @description
#'
#' Subset the `Backend` based on the provided `spectraData` data frame.
#' Subsetting could/should be done based on columns `"fileIdx"`, `"spIdx"`
#' and/or `rownames(spectraData)`.
#'
#' @param object `Backend`
#'
#' @param spectraData `DataFrame` with the spectrum metadata of the spectra to
#'     which the `object` should be subsetted.
#'
#' @return A `Backend` class.
#'
#' @author Johannes Rainer
#'
#' @rdname hidden_aliases
#'
#' @noRd
setGeneric("backendSubset", def = function(object, spectraData)
    standardGeneric("backendSubset"),
    valueClass = "Backend"
)

setMethod("backendSubset", "Backend", function(object, spectraData) {
    ## reordering of file indices is allowed
    ## (see https://github.com/lgatto/MSnbase/issues/417)
    fidx <- unique(spectraData$fileIdx)
    object@files <- object@files[fidx]
    object@modCount <- object@modCount[fidx]
    validObject(object)
    object
})

#' @description
#'
#' Split the `Backend` based on the provided `spectraData` data frame.
#' Splitting could/should be done based on columns `"fileIdx"`.
#'
#' @param object `Backend`
#'
#' @param spectraData `DataFrame` with the spectrum metadata of the spectra to
#'     which the `object` should be subsetted.
#'
#' @return A `list` of `Backend` classes.
#'
#' @author Sebastian Gibb
#'
#' @rdname hidden_aliases
#'
#' @noRd
setGeneric(
    "backendSplitByFile",
    def = function(object, spectraData, ...)
        standardGeneric("backendSplitByFile"),
    valueClass = "list"
)
setMethod("backendSplitByFile", "Backend", function(object, spectraData, ...) {
    lapply(
        split(spectraData, spectraData$fileIdx),
        backendSubset,
        object = object
    )
})

#' @description
#'
#' The reverse/undo function to `backendSplitByFile`.
#'
#' @param object `Backend`
#'
#' @param spectraData `DataFrame` with the spectrum metadata of the spectra to
#'     which the `object` should be subsetted.
#'
#' @return A `Backend` class.
#'
#' @author Sebastian Gibb
#'
#' @rdname hidden_aliases
#'
#' @noRd
setGeneric(
    "backendSplitByFile<-",
    def = function(object, spectraData, ..., value)
        standardGeneric("backendSplitByFile<-"),
    valueClass = "Backend"
)
setReplaceMethod(
    "backendSplitByFile",
    "Backend",
    function(object, spectraData, ..., value) {
    fidx <- unique(sort.int(spectraData$fileIdx))
    if (length(fidx) != length(value)) {
        stop("Length of assignment is not the same as number of files.")
    }
    for (i in seq_along(value)) {
        object@files[fidx[i]] <- value[[i]]@files
        object@modCount[fidx[i]] <- value[[i]]@modCount
    }
    validObject(object)
    object
})

#' @description
#'
#' `backendUpdateMetadata` updates the spectrum metadata on backends that
#' support it with the provided `spectraData`. It ensures that changes to the
#' metadata in the upstream object (e.g. `MSnExperiment`) are propagated to
#' the backend.
#'
#' This method is called each time the spectrum metadata is updated in the
#' `MSnExperiment`, e.g. by `spectraData(object) <- new_spd`.
#'
#' @param x `Backend`.
#'
#' @param spectraData `DataFrame` with the updated spectrum metadata.
#'
#' @return A `Backend` class.
#'
#' @author Johannes Rainer
#'
#' @rdname hidden_aliases
#'
#' @noRd
setGeneric("backendUpdateMetadata", def = function(object, spectraData)
    standardGeneric("backendUpdateMetadata"),
    valueClass = "Backend")
setMethod("backendUpdateMetadata", "Backend", function(object, spectraData) {
    object
})
