#' @include hidden_aliases.R
NULL

#' @title Mass spectrometry data managing backends
#'
#' @aliases Backend-class BackendMzR-class backendInitialize backendImportData
#'     backendReadSpectra backendWriteSpectra BackendMemory-class
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
#' experiments.
#' It mimics the classical [MSnExp-class] behaviour.
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
#'
#' It may also provide methods for:
#'
#' - [backendInitialize()] to setup the backend (create files, tables, ...).
#' - [backendImportData()] to initial import data from source *mzML* files.
#' - [backendDeepCopy()] to create a copy of the backend with associated files.
#'
#' @name Backend
#'
#' @author Sebastian Gibb
#'
#' @noRd
setClass("Backend",
    slots = c(
        files = "character"     # src files (i.e. mzML files)
    ),
    prototype = prototype(files = character()),
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
        if (is.null(names(x)))
            return("Names for 'file' missing.")
        if (anyDuplicated(names(x)))
            return("Duplicated names of 'file' found.")
        if (isFALSE(all(startsWith(names(x), "F"))))
            return("Names of 'file' don't start with 'F'.")
    }
    NULL
}

setValidity("Backend", function(object) {
    msg <- .valid.Backend.files(object@files)
    if (is.null(msg)) { TRUE } else { msg }
})

#' @rdname hidden_aliases
#' @param object Object to display.
setMethod(
    "show",
    signature="Backend",
    definition=function(object) {
    cat("Backend:", class(object)[1L], "\n")
    cat("Source files:\n",
        paste(" ", basename(object@files), collapse="\n"), "\n", sep=""
    )
})

#' @rdname hidden_aliases
setMethod("fileNames", "Backend", function(object, ...) object@files)

#' @rdname hidden_aliases
setMethod("filterFile", "Backend", function(object, file, ...) {
    if (is.character(file)) {
        file <- base::match(file, object@files)
    }
    object@files <- object@files[file]
    validObject(object)
    object
})

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
    def=function(object, files, spectraData, ...)
        standardGeneric("backendInitialize"),
    valueClass="Backend"
)

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendInitialize",
    signature="Backend",
    definition=function(object, files, spectraData, ...) {
    object@files <- normalizePath(files)
    ## use same names for files as `rownames(spectraData)`
    ## e.g. for 2 files: F1, F2; for 10 files: F01, F02, ..., F10
    names(object@files) <- sprintf(
        paste0("F%0", ceiling(log10(length(files))), "d"), seq_along(files)
    )
    validObject(object)
    object
})

#' Import spectra data into a backend
#'
#' This generic is used to import spectra data into a backend.
#'
#' It should be only reimplemented if the backend is created from scratch
#' e.g. for *HDF5* the content of the .mzML files is imported into .h5 files.
#' It must not be reimplemented for the .mzML backend.
#'
#' @inheritParams backendInitialize
#' @param BPPARAM Should parallel processing be used? See
#' [BiocParallel::bpparam()].
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setGeneric(
    "backendImportData",
    def=function(object, spectraData, ..., BPPARAM=bpparam())
        standardGeneric("backendImportData"),
    valueClass="Backend"
)

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendImportData",
    signature="Backend",
    definition=function(object, spectraData, ..., BPPARAM=bpparam()) {
    object
})

#' Create a deep copy of the backend
#'
#' This generic is used to create a deep copy of the backend and its associated
#' files.
#'
#' If an object with a file-based backend is copied in R via `newObj <- oldObj`
#' and subsequently modified the file content may be corrupted. The deep copy
#' also should copy all associated files to a new destination.
#'
#' It should be only reimplemented if the backend uses files or other resources
#' that could be corrputed when the object is copied via `<-`, e.g. *HDF5*.
#'
#' @param object A [Backend-class] derivate.
#' @inheritParams backendInitialize
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setGeneric(
    "backendDeepCopy",
    def=function(object, ...)
        standardGeneric("backendDeepCopy"),
    valueClass="Backend"
)

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendDeepCopy",
    signature="Backend",
    definition=function(object, ...) {
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
    def=function(object, spectraData, processingQueue = list(), ...)
        standardGeneric("backendReadSpectra"),
    valueClass="list"
)

#' Write spectra data to a backend
#'
#' This generic is used to write spectra data to a backend.
#'
#' It *MUST* be reimplemented by all backends!
#'
#' @inheritParams backendReadSpectra
#' @param spectra A list of [Spectrum-class] objects that should be written to
#' the backend.
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setGeneric(
    "backendWriteSpectra",
    def=function(object, spectra, spectraData, ...)
        standardGeneric("backendWriteSpectra"),
    valueClass="Backend"
)
