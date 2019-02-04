#' @include hidden_aliases.R
NULL

#' Backend class that stores data temporarily to hdf5 files.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
setClass("BackendHdf5", contains = "Backend", slots = c(md5sum = "character",
                                                        hdf5file = "character"))

setMethod("show", "BackendHdf5", function(object) {
    callNextMethod()
    cat("Hdf5 path:", unique(dirname(object@hdf5file)), "\n")
})

#' @rdname Backend
BackendHdf5 <- function() new("BackendHdf5")

.valid.BackendHdf5.md5sum <- function(x, y) {
    if (length(x) != length(y))
        "different number of md5 sums and hdf5 files"
    else NULL
}

.valid.BackendHdf5.hdf5file <- function(x, y) {
    msg <- NULL
    if (length(x) != length(y))
        msg <- "different number of hdf5 and original files"
    if (anyDuplicated(x))
        msg <- c(msg, "hdf5 files are not unique")
    if (length(msg)) msg
    else NULL
}

.valid.BackendHdf5.hdf5file.exist <- function(x) {
    if (!all(file.exists(x)))
        paste0("file(s) ", paste0(x[!file.exists(x)]), " not found")
    else NULL
}

setValidity("BackendHdf5", function(object) {
    msg <- c(.valid.BackendHdf5.md5sum(object@files, object@md5sum),
             .valid.BackendHdf5.hdf5file(object@hdf5file, object@files),
             .valid.BackendHdf5.hdf5file.exist(object@hdf5file))
    if (length(msg)) msg
    else TRUE
})

#' Initialize the Hdf5 backend, i.e. create the hdf5 files and add their names
#' to the @hdf5file slot.
#'
#' @param object `BackendHdf5`
#'
#' @param files `character` with the original file names.
#'
#' @param spectraData `DataFrame` with the spectrum metadata.
#'
#' @param path `character(1)` defining the path where the hdf5 files should be
#'     stored.
#'
#' @return `BackendHdf5` object.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
setMethod("backendInitialize", "BackendHdf5", function(object, files,
                                                       spectraData, path = ".",
                                                       ...) {
    path <- suppressWarnings(normalizePath(path))
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    n_files <- length(files)
    object@files <- files
    object@md5sum <- character(n_files)
    object@hdf5file <- character(n_files)
    comp_level <- .hdf5_compression_level()
    object@hdf5file <- file.path(path, paste0(.vdigest(files), ".h5"))
    if (any(file.exists(object@hdf5file)))
        stop("File(s) ", paste0(object@hdf5file[file.exists(object@hdf5file)],
                                collapse = ", "), "already exist(s). ",
             "Please choose a different 'path'.")
    for (i in seq_len(n_files)) {
        h5 <- H5Fcreate(object@hdf5file[i])
        h5createGroup(h5, "spectra")
        h5createGroup(h5, "md5")
        H5Fclose(h5)
    }
    validObject(object)
    object
})

#' Import the data from the raw MS files and store them to the hdf5 files.
#'
#' @inheritParams backendInitialize
#'
#' @return `BackendHdf5`
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
setMethod("backendImportData", "BackendHdf5", function(object, spectraData,
                                                       ...,
                                                       BPPARAM = bpparam()) {
    object@md5sum <- bpmapply(fileNames(object), object@hdf5file,
                              FUN = .serialize_msfile_to_hdf5,
                              BPPARAM = BPPARAM)
    validObject(object)
    object
})

## setMethod("backendReadSpectra", "BackendHdf5", function(object, spectraData,
##                                                         ...) {
## })

## setMethod("backendWriteSpectra", "BackendHdf5", function(object, spectra,
##                                                          spectraData) {
## })

#' Write the content of a single mzML/etc file to an h5file. We're using the
#' spectrum index in the file as data set ID.
#'
#' @return `character(1)` with the md5 sum of the spectra data.
#'
#' @author Johannes Rainer
#'
#' @noRd
.serialize_msfile_to_hdf5 <- function(file, h5file) {
    h5 <- H5Fopen(h5file)
    on.exit(H5Fclose(h5))
    comp_level <- .hdf5_compression_level()
    fh <- openMSfile(file)
    hdr <- header(fh)
    pks <- peaks(fh)
    close(fh)
    for (i in seq_along(pks)) {
        .pks <- pks[[i]]
        colnames(.pks) <- c("mz", "intensity")
        h5write(.pks, h5, paste0("spectra/", i), level = comp_level)
        pks[[i]] <- .pks
    }
    pks_md5 <- digest(pks)
    h5write(pks_md5, h5, paste0("md5/md5"), level = comp_level)
    pks_md5
}
