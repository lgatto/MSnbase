#' @include hidden_aliases.R
NULL

#' Backend class that stores data temporarily to hdf5 files.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
setClass("BackendHdf5", contains = "Backend", slots = c(checksums = "character",
                                                        h5files = "character"))

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

## setMethod("backendReadSpectra", "BackendHdf5", function(object, spectraData,
##                                                         ...) {
## })

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
#' @param fdata `data.frame` representing the feature data (`fData`) of the
#'     spectra to be returned.
#'
#' @param h5file `character`: the name of the HDF5 file.
#'
#' @param h5md5 `character`: the md5 sum of the spectrum data in a hdf5 file.
#'     This is checked against the hdf5 from the file and an error if thrown
#'     if they differ, i.e. if the data was changed in the HDF5 file.
#'
#' @param file_name `character` with the file name of the original mzML file.
#'
#' @return list of `Spectrum` objects in the order of the spectra given in
#'     param `fdata`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_read_spectra <- function(fdata, h5file, h5md5, file_name) {
    suppressPackageStartupMessages(require(MSnbase, quietly = TRUE))
    fid <-.Call("_H5Fopen", h5file, 0L, PACKAGE = "rhdf5")
    on.exit(invisible(.Call("_H5Fclose", fid, PACKAGE = "rhdf5")))
    grp_name <- paste0(.h5_group_name(file_name), "/")
    pks_md5 <- .h5_read_bare(fid, paste0(grp_name, "md5"))
    if (pks_md5 != h5md5)
        stop("The data in the Hdf5 files associated with this object appear ",
             "to have changed! Please see the Notes section in ?Hdf5MSnExp ",
             "for more information.")
    mzi <- base::lapply(paste0(grp_name, fdata$spIdx),
                        .h5_read_bare, file = fid)
    res <- vector("list", nrow(fdata))
    names(res) <- rownames(fdata)
    ms1 <- which(fdata$msLevel == 1)
    n_peaks <- base::lengths(mzi, use.names = FALSE) / 2
    if (length(ms1)) {
        mzi_ms1 <- do.call(rbind, mzi[ms1])
        res[ms1] <- Spectra1_mz_sorted(
            peaksCount = n_peaks[ms1],
            rt = fdata$retentionTime[ms1],
            acquisitionNum = fdata$acquisitionNum[ms1],
            scanIndex = fdata$spIdx[ms1],
            tic = fdata$totIonCurrent[ms1],
            mz = mzi_ms1[, 1],
            intensity = mzi_ms1[, 2],
            fromFile = fdata$fileIdx[ms1],
            centroided = fdata$centroided[ms1],
            smoothed = fdata$smoothed[ms1],
            polarity = fdata$polarity[ms1],
            nvalues = n_peaks[ms1])
    }
    msn <- which(fdata$msLevel > 1)
    if (length(msn)) {
        mzi_msn <- do.call(rbind, mzi[msn])
        res[msn] <- Spectra2_mz_sorted(
            msLevel = fdata$msLevel[msn],
            peaksCount = n_peaks[msn],
            rt = fdata$retentionTime[msn],
            acquisitionNum = fdata$acquisitionNum[msn],
            scanIndex = fdata$spIdx[msn],
            tic = fdata$totIonCurrent[msn],
            mz = mzi_msn[, 1],
            intensity = mzi_msn[, 2],
            fromFile = fdata$fileIdx[msn],
            centroided = fdata$centroided[msn],
            smoothed = fdata$smoothed[msn],
            polarity = fdata$polarity[msn],
            merged = fdata$mergedScan[msn],
            precScanNum = fdata$precursorScanNum[msn],
            precursorMz = fdata$precursorMZ[msn],
            precursorIntensity = fdata$precursorIntensity[msn],
            precursorCharge = fdata$precursorCharge[msn],
            collisionEnergy = fdata$collisionEnergy[msn],
            nvalues = n_peaks[msn])
    }
    res
}

## setMethod("backendWriteSpectra", "BackendHdf5", function(object, spectra,
##                                                          spectraData) {
## })
