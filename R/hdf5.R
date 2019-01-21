validHdf5MSnExp <- function(object, check_md5 = TRUE) {
    msg <- character()
    if (length(object)) {
        if (length(object@hdf5file) != length(fileNames(object)))
            msg <- c(msg, "Number of hdf5 files and raw files does not match")
        if (!all(file.exists(object@hdf5file)))
            msg <- c(
                msg, paste0("Hdf5 files ",
                            paste(object@hdf5file[!file.exists(object@hdf5file)],
                                  collapse = ", "), " do not exist"))
        if (length(object@hdf5file) != length(object@md5sum))
            msg <- c(msg, "Not for all Hdf5 files an MD5 sum is abailable.")
        if (check_md5 && !.h5_check_md5(object))
            msg <- c(msg, "MD5 sum in the object and the hdf5 file(s) differ(s).")
    }
    if (length(msg))
        msg
    else TRUE
}

#' Compare the md5 sum stored in the object with the sum in the hdf5 file(s).
#'
#' @param object `Hdf5MSnExp`
#'
#' @return `TRUE` if md5 sums match for all files or `FALSE` otherwise.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @md
.h5_check_md5 <- function(object) {
    res <- mapply(fileNames(object), object@hdf5file, object@md5sum,
                  FUN = function(fl, h5fl, md) {
                      h5 <- H5Fopen(h5fl)
                      md_fl <- h5read(h5, paste0(.h5_group_name(fl), "/md5"))
                      H5Fclose(h5)
                      md_fl == md
                  })
    all(unlist(res))
}

#' @title The `Hdf5MSnExp` Class for MS Data And Meta-Data
#'
#' @description
#'
#' The `Hdf5MSnExp` class encapsulates data and meta-data for mass
#' spectrometry experiments like the `MSnExp` and `OnDiskMSnExp`
#' classes. `Hdf5MSnExp` implements, like `OnDiskMSnExp`, an
#' *on disk* model, where the raw data (M/Z and intensities) are
#' stored on disk (rather than in memory, like `MSnExp` data) using,
#' as the name implies, the hdf5 file format. See [`MSnExp`] and
#' [`OnDiskMSnExp`] for a more general description of the data,
#' slots and operations.
#'
#' Object of the class are currently created with the
#' `readHdf5DiskMSData` function. Later, this backend will be
#' included in the main `readMSData` function.
#'
#' @slot hdf5file `character` containing the hdf5 filenames, one file
#'     per input mzML file/sample.
#'
#' @rdname Hdf5MSnExp-class
#'
#' @aliases Hdf5MSnExp
#'
#' @seealso See [`MSnExp`] and [`OnDiskMSnExp`] classes.
#'
#' @author Laurent Gatto, Johannes Rainer
#'
#' @md
#'
#' @examples
#' f <- msdata::proteomics(pattern = "MS3TMT11", full.names = TRUE)
#' x <- readHdf5DiskMSData(f, hdf5path = tempdir())
#' x
#'
#' ## automatically generated hdf5 file
#' hdf5FileName(x)
#'
#' x[[1]]
#' x[[2]]
#' x[[10]]
#' filterMsLevel(x, 3L)[[1]]
#'
#' ## Perform data manipulations, remove all peaks with an intensity
#' ## smaller than 1000
#' x_clean <- clean(removePeaks(x, t = 1000), all = TRUE)
#'
#' x_clean[[12]]
#'
#' ## Data manipulations are, as for OnDiskMSnExp objects, performed
#' ## on-the-fly. The data in the hdf5 files was not changed by the
#' ## clean and removePeaks calls above. Object x is thus still valid
#' ## and we can still use it.
#' MSnbase:::validHdf5MSnExp(x)
#' x[[12]]
#'
#' ## updateHdf5Data will make all data manipulations persistent, it will
#' ## overwrite the data in the hdf5 files associated with x_clean and
#' ## accessing the data afterwards will no longer require to apply
#' ## data manipulations on-the-fly.
#' updateHdf5Data(x_clean)
#'
#' x_clean[[12]]
#'
#' ## Note that the original object x is however still associated with the
#' ## same hdf5 files. This object will now no longer be valid, as the data
#' ## has changed for it.
#' MSnbase:::validHdf5MSnExp(x)
#'
#' ## Accessing the data such as with x[[12]] will thus result in an error.
#' ## The data can however be *restored* with the `convertToHdf5MSnExp`
#' ## function that will re-read all data into hdf5 files. Note that a
#' ## different `hdf5path` has to be used.
#' x <- convertToHdf5MSnExp(x, hdf5path = paste0(tempdir(), "/2/"))
#'
#' MSnbase:::validHdf5MSnExp(x)
#' x[[12]]
#'
#' ## clean up session
#' file.remove(hdf5FileName(x))
.Hdf5MSnExp <- setClass("Hdf5MSnExp",
                        slots = c(hdf5file = "character",
                                  md5sum = "character"),
                        contains = "OnDiskMSnExp",
                        prototype = prototype(
                            new("VersionedBiobase",
                                versions = c(classVersion("MSnExp"),
                                             classVersion("OnDiskMSnExp"),
                                             Hdf5MSnExp = "0.0.1")),
                            spectraProcessingQueue = list(),
                            backend = character()),
                        validity = function(object) validHdf5MSnExp(object, FALSE))

## .onDisk2hdf5 <- function(from, filename) {
##           to <- .Hdf5MSnExp()
##           for (sl in .slotNames0(from))
##               slot(to, sl) <- slot(from, sl)
##           to@hdf5file <- filename
##           to
##       }

#' Write the content of a single mzML/etc file to an h5file. We're using the
#' spectrum index in the file as data set ID.
#'
#' @return `character(1)` with the md5 sum of the spectra data.
#'
#' @author Johannes Rainer
#'
#' @noRd
.serialize_msfile_to_hdf5 <- function(file, h5file) {
    h5 <- H5Fcreate(h5file)
    comp_level <- .hdf5_compression_level()
    fh <- openMSfile(file)
    hdr <- header(fh)
    pks <- peaks(fh)
    grp_name <- .h5_group_name(file)
    h5createGroup(h5, grp_name)
    grp_name <- paste0(grp_name, "/")
    close(fh)
    for (i in seq_along(pks)) {
        .pks <- pks[[i]]
        colnames(.pks) <- c("mz", "intensity")
        h5write(.pks, h5, paste0(grp_name, i), level = comp_level)
        pks[[i]] <- .pks
    }
    pks_md5 <- digest::digest(pks)
    h5write(pks_md5, h5, paste0(grp_name, "md5"), level = comp_level)
    H5Fclose(h5)
    pks_md5
}

serialise_to_hdf5 <- function(object, path, BPPARAM = bpparam()) {
    stopifnot(inherits(object, "OnDiskMSnExp"))
    path <- suppressWarnings(normalizePath(path))
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    filename <- paste0(path, "/", sapply(fileNames(object), digest::sha1), ".h5")
    if (any(file.exists(filename)))
        stop("File(s) ", paste(filename[file.exists(filename)], collapse = ", "),
             " already exist(s). Please use a different 'hdf5path'.")
    md5s <- bpmapply(fileNames(object), filename,
                     FUN = .serialize_msfile_to_hdf5,
                     BPPARAM = BPPARAM)
    object@hdf5file <- filename
    object@md5sum <- md5s
    object
}

#' @rdname Hdf5MSnExp-class
#'
#' @description
#'
#' `convertToHdf5MSnExp` coerces an [OnDiskMSnExp] object (or any object
#' extending it) to an `Hdf5MSnExp` object.
#'
#' @md
convertToHdf5MSnExp <- function(object, verbose = FALSE, hdf5path = ".",
                                BPPARAM = bpparam()) {
    stopifnot(inherits(object, "OnDiskMSnExp"))
    object <- as(object, "Hdf5MSnExp")
    if (verbose) message("Serialising to hdf5...")
    object <- serialise_to_hdf5(object, hdf5path, BPPARAM = BPPARAM)
    if (validObject(object)) object
}

#' @rdname Hdf5MSnExp-class
readHdf5DiskMSData <- function(files, pdata = NULL, msLevel. = NULL,
                               verbose = isMSnbaseVerbose(),
                               centroided. = NA, smoothed. = NA,
                               hdf5path = ".", BPPARAM = bpparam()) {
    obj <- readOnDiskMSData(normalizePath(files), pdata, msLevel.,
                            verbose, centroided.,
                            smoothed.)
    convertToHdf5MSnExp(obj, verbose = verbose, hdf5path = hdf5path,
                        BPPARAM = BPPARAM)
}

## ##' The `hdf5Close` and `hdf5Open` function respectively close and
## ##' open the connection to the `Hdf5MSnExp` object hdf5 file. The both
## ##' return the input object with the updated hdf5 file
## ##' handle. `isHdf5Open` test if the hdf5 file is open and retuns a
## ##' `logical(1)`.
## ##'
## ##' @param object An instance of class `Hdf5MSnExp`.
## ##' @rdname Hdf5MSnExp-class
## hdf5Close <- function(object) {
##     if (isHdf5Open(object))
##         tryCatch(rhdf5::H5Fclose(object@hdf5handle),
##                  error = function(e)
##                      stop("Encountered error trying to close ",
##                           object@hdf5file,
##                           ". Use `h5closeAll()` to force close all hdf5 files."))
##     invisible(TRUE)
## }

## ##' @rdname Hdf5MSnExp-class
## hdf5Open <- function(object) {
##     if (!isHdf5Open(object))
##         object@hdf5handle <- rhdf5::H5Fopen(object@hdf5file)
##     object
## }

## ##' @rdname Hdf5MSnExp-class
## isHdf5Open <- function(object) {
##     stopifnot(inherits(object, "Hdf5MSnExp"))
##     if (is.null(object@hdf5handle))
##         return(FALSE)
##     rhdf5::H5Iis_valid(object@hdf5handle)
## }

#' @rdname Hdf5MSnExp-class
hdf5FileName <- function(object) {
    stopifnot(inherits(object, "Hdf5MSnExp"))
    object@hdf5file
}

## ##' @rdname Hdf5MSnExp-class
## `hdf5FileName<-` <- function(object, value) {
##     stopifnot(inherits(object, "Hdf5MSnExp"))
##     stopifnot(inherits(value, "character"))
##     if (length(value) != 1)
##         stop("Please provide a single file name")
##     if (isHdf5Open(object))
##         object <- hdf5Close(object)
##     stopifnot(file.rename(hdf5FileName(object), value))
##     object@hdf5file <- value
##     msg <- validHdf5MSnExp(object)
##     if (is.character(msg)) stop(msg)
##     object
## }

.h5_group_name <- function(x) {
    vapply(x, digest::sha1, character(1), USE.NAMES = FALSE)
}

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
#' @param h5file `character(1)` with the name of the hdf5 file into which the
#'     data should be written.
#'
#' @param group `character(1)` with the hdf5 group name.
#'
#' @param prune `logical(1)` whether the group should be deleted before writing
#'     the data (see description above for more details).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.h5_write_spectra <- function(x, h5file, group, prune = TRUE) {
    h5 <- H5Fopen(h5file)
    on.exit(invisible(H5Fclose(h5)))
    comp_level <- .hdf5_compression_level()
    if (prune) {
        h5delete(h5, group)
    }
    if (!H5Lexists(h5, group))
        h5createGroup(h5, group)
    for (i in seq_along(x)) {
        h5write(cbind(mz = mz(x[[i]]), intensity = intensity(x[[i]])),
                h5, name = paste0(group, "/", names(x)[i]),
                level = comp_level)
    }
    ## Create MD5 on spectra (instead of list of matrices) and write
    h5_md5 <- digest(x)
    h5write(h5_md5, file = h5, name = paste0(group, "/md5"), level = comp_level)
    h5_md5
}

#' @rdname Hdf5MSnExp-class
setMethod("spectrapply", "Hdf5MSnExp", function(object, FUN = NULL,
                                                BPPARAM = bpparam(), ...) {
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)
    isOK <- validateFeatureDataForOnDiskMSnExp(fData(object))
    if (!is.null(isOK))
        stop(isOK)
    fDataPerFile <- split.data.frame(fData(object),
                                     f = fData(object)$fileIdx)
    h5_fls <- object@hdf5file
    h5_md5 <- object@md5sum
    theQ <- processingQueue(object)
    if (!is.null(FUN))
        theQ <- c(theQ, list(ProcessingStep(FUN, ARGS = list(...))))
    vals <- bpmapply(fDataPerFile, h5_fls, h5_md5, fileNames(object),
                     FUN = function(fdata, h5_file, md5, file_name, queue) {
                         .apply_processing_queue(
                             .h5_read_spectra(fdata, h5_file, md5, file_name),
                             queue)
                     }, MoreArgs = list(queue = theQ),
                     BPPARAM = BPPARAM, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    vals <- unlist(vals, recursive = FALSE)
    vals[rownames(fData(object))]
})

setMethod("filterFile", "Hdf5MSnExp", function(object, file) {
    if (missing(file))
        return(object)
    if (is.character(file)) {
        file <- base::match(file, basename(fileNames(object)))
    }
    file <- base::sort(unique(file))
    object@hdf5file <- object@hdf5file[file]
    object@md5sum <- object@md5sum[file]
    callNextMethod(object, file)
})

#' @description
#'
#' `updateHdf5Data` *consolidates* an `Hdf5MSnExp` object by applying all
#' registered processing steps to each spectrum and saving the updated
#' data to the hdf5 file(s) associated with the sample(s) in `x`. This has
#' some implications on copies of the input object. See notes for more
#' information.
#'
#' @note
#'
#' `updateHdf5Data` will overwrite the content of the hdf5 file(s) associated
#' with the `Hdf5MSnExp`. Thus, if another copy of the object, prior to any data
#' manipulations, exists that points to the same hdf5 files, that object might
#' become corrupt. Note that `convertToHdf5MSnExp` can *restore* corrupted
#' `Hdf5MSnExp` files again (see examples for details).
#'
#' @param x `Hdf5MSnExp` object.
#'
#' @md
#'
#' @rdname Hdf5MSnExp-class
updateHdf5Data <- function(x) {
    if (isMSnbaseVerbose())
        message("Note: after writing data to hdf5 files any copy of the input",
                " object linking to the same hdf5 files will no longer be ",
                "valid.\n")
    stopifnot(inherits(x, "Hdf5MSnExp"))
    eval.parent(call("<-", substitute(x), .consolidate(x)), n = 1)
}
.consolidate <- function(x, BPPARAM = bpparam(), ...) {
    x_split <- splitByFile(x, f = factor(fileNames(x)))
    if (isMSnbaseVerbose())
        message("Consolidating ", length(x_split), " files.")
    md5s <- bplapply(x_split, function(z) {
        suppressPackageStartupMessages(require(MSnbase, quietly = TRUE))
        grp <- .h5_group_name(fileNames(z))
        sps <- spectra(z)
        names(sps) <- fData(z)$spIdx
        .h5_write_spectra(sps, z@hdf5file, group = grp)
    }, BPPARAM = BPPARAM)
    x@spectraProcessingQueue <- list()
    x@md5sum <- unlist(md5s)
    x
}
