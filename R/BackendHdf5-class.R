#' @include hidden_aliases.R
NULL

#' Backend class that stores data temporarily to hdf5 files.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @md
#'
#' @noRd
setClass("BackendHdf5", contains = "Backend", slots = c(checksums = "character",
                                                        h5files = "character"))

#' @rdname hidden_aliases
setMethod("show", "BackendHdf5", function(object) {
    callNextMethod()
    cat("Hdf5 path:", unique(dirname(object@h5files)), "\n")
})

#' @rdname Backend
BackendHdf5 <- function() new("BackendHdf5")

.valid.BackendHdf5.checksums <- function(x, y) {
    if (length(x) != length(y))
        "different number of md5 sums and hdf5 files"
    else NULL
}

.valid.BackendHdf5.h5files <- function(x, y) {
    msg <- NULL
    if (length(x) != length(y))
        msg <- "different number of hdf5 and original files"
    if (anyDuplicated(x))
        msg <- c(msg, "hdf5 files are not unique")
    if (length(msg)) msg
    else NULL
}

.valid.BackendHdf5.h5files.exist <- function(x) {
    if (!all(file.exists(x)))
        paste0("file(s) ", paste0(x[!file.exists(x)]), " not found")
    else NULL
}

setValidity("BackendHdf5", function(object) {
    msg <- c(.valid.BackendHdf5.checksums(object@files, object@checksums),
             .valid.BackendHdf5.h5files(object@h5files, object@files),
             .valid.BackendHdf5.h5files.exist(object@h5files))
    if (length(msg)) msg
    else TRUE
})

#' Initialize the Hdf5 backend, i.e. create the hdf5 files and add their names
#' to the @h5files slot.
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
    names(object@files) <- sprintf(
        paste0("F%0", ceiling(log10(length(files))), "d"), seq_along(files)
    )
    object@checksums <- character(n_files)
    object@h5files <- character(n_files)
    comp_level <- .hdf5_compression_level()
    object@h5files <- file.path(path, paste0(.vdigest(files), ".h5"))
    if (any(file.exists(object@h5files)))
        stop("File(s) ", paste0(object@h5files[file.exists(object@h5files)],
                                collapse = ", "), "already exist(s). ",
             "Please choose a different 'path'.")
    for (i in seq_len(n_files)) {
        h5 <- H5Fcreate(object@h5files[i])
        h5createGroup(h5, "spectra")
        h5createGroup(h5, "checksum")
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
    object@checksums <- bpmapply(fileNames(object), object@h5files,
                                 FUN = .serialize_msfile_to_hdf5,
                                 BPPARAM = BPPARAM)
    validObject(object)
    object
})

#' Write the content of a single mzML/etc file to an h5file. We're using the
#' spectrum index in the file as data set ID.
#'
#' @return `character(1)` with the md5 sum of the spectra data.
#'
#' @author Johannes Rainer, Sebastian Gibb
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
    spids <- paste0("/spectra/", seq_along(pks))
    for (i in seq_along(pks)) {
        h5write(pks[[i]], h5, spids[i], level = comp_level)
    }
    H5Fflush(h5)
    checksum <- digest(h5file, file = TRUE)
    h5write(checksum, h5, paste0("/checksum/checksum"), level = comp_level)
    checksum
}

setMethod("backendReadSpectra", "BackendHdf5", function(object, spectraData,
                                                        ...) {
    file_f <- factor(spectraData$fileIdx, levels = unique(spectraData$fileIdx))
    idx <- as.integer(levels(file_f))
    fls <- object@h5files[idx]
    checksums <- object@checksums[idx]
    if (length(fls) == 1)
        .h5_read_spectra(spectraData, fls, checksums)
    else
        unlist(mapply(split(spectraData, file_f), fls, checksums,
                      FUN = .h5_read_spectra, SIMPLIFY = FALSE,
                      USE.NAMES = FALSE), recursive = FALSE)
})

#' This is based on Mike Smith's function (issue #395) that directly accesses
#' the data without validation and file checking.
#'
#' @param file the ID as returned by `_H5Fopen`.
#'
#' @param name `character` defining the data set to be read.
#'
#' @return the imported data set (in most cases a `matrix`).
#'
#' @author Mike Smith, Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_read_bare <- function(file, name = "") {
    did <- .Call("_H5Dopen", file, name, NULL, PACKAGE = "rhdf5")
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, FALSE,
                 PACKAGE = "rhdf5")
    invisible(.Call("_H5Dclose", did, PACKAGE = "rhdf5"))
    res
}

#' Read spectrum data from an hdf5 file and return a list of `Spectrum` objects.
#'
#' @note
#'
#' This function uses a constructor function that creates all `Spectrum`
#' objects in C++ for added performance.
#'
#' @param spectraData `DataFrame` representing the spectrum metadata (for
#'     spectra from a single file).
#'
#' @param h5file `character`: the name of the HDF5 file.
#'
#' @param checksum `character`: the checksum of the hdf5 file.
#'     This is checked against the hdf5 from the file and an error if thrown
#'     if they differ, i.e. if the data was changed in the HDF5 file.
#'
#' @return list of `Spectrum` objects in the order of the spectra given in
#'     param `fdata`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_read_spectra <- function(spectraData, h5file, checksum) {
    suppressPackageStartupMessages(require(MSnbase, quietly = TRUE))
    fid <- .Call("_H5Fopen", h5file, 0L, PACKAGE = "rhdf5")
    on.exit(invisible(.Call("_H5Fclose", fid, PACKAGE = "rhdf5")))
    h5_checksum <- .h5_read_bare(fid, "/checksum/checksum")
    if (h5_checksum != checksum)
        stop("The data in the Hdf5 files associated with this object appear ",
             "to have changed! Please see the Notes section in ?Backend ",
             "for more information.")
    .spectra_from_data(base::lapply(paste0("/spectra/", spectraData$spIdx),
                                    .h5_read_bare, file = fid), spectraData)
}

#' Write spectra of an mzML file to a group within a hdf5 file. Depending on
#' the value of `prune` the hdf5 group will be deleted before writing the data.
#' This is required as datasets in a hdf5 file can not be *updated* if their
#' dimensions don't match.
#'
#' @note
#'
#' `Hdf5MSnExp` object store spectrum data in (custom) hdf5 files for faster
#' data access with data for each sample (from each input file) being stored in
#' its own hdf5 file. Copies of `Hdf5MSnExp` objects, such as created with
#' `a <- b` where `b` is an `Hdf5MSnExp` object will be associated with the
#' same set of hdf5 files.
#'
#' Each spectrum will be saved as a dataset with the name being the name of the
#' Spectrum in `x`. The names of `x` should thus be the **spectrum index** in
#' the original file (i.e. `fData(object)$spIdx`) and not the spectrum names
#' such as `"F01.001"`!
#'
#' @param x `list` of `Spectrum` objects.
#'
#' @param spectraData `DataFrame` with the spectra metadata information.
#'
#' @param h5file `character(1)` with the name of the hdf5 file into which the
#'     data should be written.
#'
#' @param prune `logical(1)` whether the group should be deleted before writing
#'     the data (see description above for more details).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_write_spectra <- function(x, spectraData, h5file, prune = TRUE) {
    h5 <- H5Fopen(h5file)
    on.exit(invisible(H5Fclose(h5)))
    comp_level <- .hdf5_compression_level()
    if (prune) {
        h5delete(h5, "spectra")
    }
    if (!H5Lexists(h5, "spectra"))
        h5createGroup(h5, "spectra")
    spids <- paste0("/spectra/", spectraData$spIdx)
    for (i in seq_along(x)) {
        h5write(cbind(mz = mz(x[[i]]), intensity = intensity(x[[i]])),
                h5, name = spids[i], level = comp_level)
    }
    H5Fflush(h5)
    checksum <- digest(h5file, file = TRUE)
    h5write(checksum, h5, "/checksum/checksum", level = comp_level)
    checksum
}

setMethod("backendWriteSpectra", "BackendHdf5", function(object, spectra,
                                                         spectraData) {
    file_f <- factor(spectraData$fileIdx, levels = unique(spectraData$fileIdx))
    idx <- as.integer(levels(file_f))
    fls <- object@h5files[idx]
    if (length(fls) == 1)
        chk <- .h5_write_spectra(spectra, spectraData, fls)
    else
        chk <- unlist(mapply(split(spectra, file_f), split(spectraData, file_f),
                             fls, FUN = .h5_write_spectra,
                             SIMPLIFY = FALSE, USE.NAMES = FALSE),
                      recursive = FALSE)
    object@checksums[idx] <- chk
    object
})

setMethod("backendSubset", signature(object = "BackendHdf5", i = "numeric",
                                     file = "numeric"),
          function(object, i, file) {
              object@files <- object@files[file]
              object@checksums <- object@checksums[file]
              object@h5files <- object@h5files[file]
              validObject(object)
              object
          })
