#' @include hidden_aliases.R
NULL

#' Backend class that stores data temporarily to hdf5 files.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @md
#'
#' @noRd
setClass("BackendHdf5",
    contains = "Backend",
    slots = c(
        h5files = "character"
    )
)

#' @rdname hidden_aliases
setMethod("show", "BackendHdf5", function(object) {
    callNextMethod()
    cat("Hdf5 path:", unique(dirname(object@h5files)), "\n")
})

#' @rdname Backend
BackendHdf5 <- function() new("BackendHdf5")

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
    msg <- c(.valid.BackendHdf5.h5files(object@h5files, object@files),
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
    path <- normalizePath(path, mustWork = FALSE)
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    comp_level <- .hdf5_compression_level()
    object@h5files <- file.path(path, paste0(.vdigest(files), ".h5"))
    if (any(file.exists(object@h5files)))
        stop("File(s) ", paste0(object@h5files[file.exists(object@h5files)],
                                collapse = ", "), " already exist(s). ",
             "Please choose a different 'path'.")
    for (i in seq_along(object@h5files)) {
        h5 <- H5Fcreate(object@h5files[i])
        h5createGroup(h5, "spectra")
        h5createGroup(h5, "modification")
        H5Fclose(h5)
    }
    callNextMethod()
})

setMethod("backendReadSpectra", "BackendHdf5", function(object, spectraData,
                                                        ...) {
    file_f <- factor(spectraData$fileIdx, levels = unique(spectraData$fileIdx))
    idx <- as.integer(levels(file_f))
    fls <- object@h5files[idx]
    modCount <- object@modCount[idx]
    if (length(fls) == 1)
        .h5_read_spectra(spectraData, fls, modCount)
    else
        unlist(mapply(split(spectraData, file_f), fls, modCount,
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
#' @param modCount `integer`: the modification counter stored in the R object
#'     This is checked against one stored in the hdf5 file and an error is thrown
#'     if they differ, i.e. if the object was copied and the hdf5 files modified
#'     by the original/source object.
#'
#' @return list of `Spectrum` objects in the order of the spectra given in
#'     param `fdata`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_read_spectra <- function(spectraData, h5file, modCount) {
    suppressPackageStartupMessages(require(MSnbase, quietly = TRUE))
    fid <- .Call("_H5Fopen", h5file, 0L, PACKAGE = "rhdf5")
    on.exit(invisible(.Call("_H5Fclose", fid, PACKAGE = "rhdf5")))
    h5modCount <- .h5_read_bare(fid, "/modification/counter")
    if (h5modCount != modCount)
        stop("The data in the hdf5 files associated with this object appear ",
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
#' `BackendHdf5` object store spectrum data in (custom) hdf5 files for faster
#' data access with data for each sample (from each input file) being stored in
#' its own hdf5 file. Copies of `BackendHdf5` objects, such as created with
#' `a <- b` where `b` is an `BackendHdf5` object will be associated with the
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
#' @param modCount `integer`, modification counter, incremented by one for each
#'     write
#'
#' @param prune `logical(1)` whether the group should be deleted before writing
#'     the data (see description above for more details).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_write_spectra <- function(x, spectraData, h5file, modCount, prune = TRUE) {
    h5 <- H5Fopen(h5file)
    on.exit(invisible(H5Fclose(h5)))
    comp_level <- .hdf5_compression_level()
    x <- force(x)
    if (prune)
        h5delete(h5, "/spectra")
    if (!H5Lexists(h5, "/spectra"))
        h5createGroup(h5, "/spectra")
    spids <- paste0("/spectra/", spectraData$spIdx)
    for (i in seq_along(x)) {
        h5write(cbind(mz(x[[i]]), intensity(x[[i]])),
                h5, name = spids[i], level = comp_level)
    }
    h5write(modCount, h5, "/modification/counter", level = comp_level)
}

setMethod("backendWriteSpectra", "BackendHdf5", function(object, spectra,
                                                         spectraData,
                                                         updateModCount, ...) {
    file_f <- factor(spectraData$fileIdx, levels = unique(spectraData$fileIdx))
    idx <- as.integer(levels(file_f))
    fls <- object@h5files[idx]
    if (updateModCount)
        object@modCount[idx] <- object@modCount[idx] + 1L
    if (length(fls) == 1)
        .h5_write_spectra(spectra, spectraData, fls, object@modCount[idx])
    else
        unlist(mapply(split(spectra, file_f), split(spectraData, file_f),
                            fls, object@modCount[idx], FUN = .h5_write_spectra,
                            SIMPLIFY = FALSE, USE.NAMES = FALSE),
                      recursive = FALSE)
    object
})

setMethod("backendSubset", "BackendHdf5", function(object, spectraData) {
    object@h5files <- object@h5files[unique(spectraData$fileIdx)]
    callNextMethod()
})

setReplaceMethod(
    "backendSplitByFile",
    "BackendHdf5",
    function(object, spectraData, ..., value) {
    fidx <- unique(sort.int(spectraData$fileIdx))
    for (i in seq(along=value)) {
        object@h5files[fidx[i]] <- value[[i]]@h5files
    }
    callNextMethod()
})
