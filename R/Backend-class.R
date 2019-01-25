#' @include hidden_aliases.R
NULL

#' Backend class
#'
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
#' @name Backend-class
#' @docType class
#' @family Backend classes
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @export
setClass("Backend", contains="VIRTUAL")

#' @rdname hidden_aliases
#' @param object Object to display.
#' @export
setMethod(
    "show",
    signature="Backend",
    definition=function(object) {
    cat("Backend:", class(object)[1L], "\n")
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
#' @param ... Other arguments passed to the methods.
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @export
setGeneric(
    "backendInitialize",
    def=function(object, files, ...) standardGeneric("backendInitialize"),
    valueClass="Backend"
)

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendInitialize",
    signature="Backend",
    definition=function(object, files, ...) {
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
#' @param spectraData A [S4Vectors::DataFrame-class]
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @export
setGeneric(
    "backendImportData",
    def=function(object, files, spectraData, ..., BPPARAM=bpparam())
        standardGeneric("backendImportData"),
    valueClass="Backend"
)

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendImportData",
    signature="Backend",
    definition=function(object, files, spectraData, ..., BPPARAM=bpparam()) {
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
#' @inheritParams backendImportData
#' @return A [Backend-class] derivate.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @export
setGeneric(
    "backendDeepCopy",
    def=function(object, ..., BPPARAM=bpparam())
        standardGeneric("backendDeepCopy"),
    valueClass="Backend"
)

#' @rdname hidden_aliases
#' @export
setMethod(
    "backendDeepCopy",
    signature="Backend",
    definition=function(object, ..., BPPARAM=bpparam()) {
    object
})

#' Read spectra data from a backend
#'
#' This generic is used to read spectra data from a backend.
#'
#' It *MUST* be reimplemented by all backends!
#'
#' @inheritParams backendImportData
#' @param file The path to the source (generally .mzML) file.
#' @return A list of [Spectrum-class] objects.
#' @family Backend generics
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @export
setGeneric(
    "backendReadSpectra",
    def=function(object, file, spectraData, ..., BPPARAM=bpparam())
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
#' @export
setGeneric(
    "backendWriteSpectra",
    def=function(object, file, spectra, spectraData, ..., BPPARAM=bpparam())
        standardGeneric("backendWriteSpectra"),
    valueClass="Backend"
)
