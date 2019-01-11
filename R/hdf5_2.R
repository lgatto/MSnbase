.Hdf5MSnExp2 <- setClass("Hdf5MSnExp2",
                        slots = c(hdf5file = "character"),
                        contains = "OnDiskMSnExp",
                        prototype = prototype(
                            new("VersionedBiobase",
                                versions = c(classVersion("MSnExp"),
                                             classVersion("OnDiskMSnExp"),
                                             Hdf5MSnExp2 = "0.0.1")),
                            spectraProcessingQueue = list(),
                            backend = character()))

.onDisk2hdf52 <- function(from, filename) {
          to <- .Hdf5MSnExp2()
          for (sl in .slotNames0(from))
              slot(to, sl) <- slot(from, sl)
          to@hdf5file <- filename
          to
      }

validHdf5MSnExp2 <- function(object) {
    msg <- character()
    if (length(object)) {
        if (!length(object@hdf5file) == length(fileNames(object)))
            msg <- c(msg, "Number of hdf5 files and raw files does not match")
        if (!all(file.exists(object@hdf5file)))
            msg <- c(
                msg, paste0("Hdf5 files ",
                            paste(object@hdf5file[!file.exists(object@hdf5file)],
                                  collapse = ", "), " do not exist"))
    }
    if (length(msg))
        msg
    else TRUE
}

#' Write the content of a single mzML/etc file to an h5file. We're using the
#' spectrum index in the file as data set ID.
#'
#' @noRd
.serialize_msfile_to_hdf5 <- function(file, h5file) {
    h5 <- rhdf5::H5Fcreate(h5file)
    comp_level <- .hdf5_compression_level()
    fh <- openMSfile(file)
    hdr <- header(fh)
    pks <- peaks(fh)
    close(fh)
    for (i in seq_along(pks)) {
        .pks <- pks[[i]]
        colnames(.pks) <- c("mz", "intensity")
        rhdf5::h5write(.pks, h5, as.character(i), level = comp_level)
    }
    rhdf5::H5Fclose(h5)
    invisible(h5file)
}

serialise_to_hdf52 <- function(object, path, BPPARAM = bpparam()) {
    stopifnot(inherits(object, "OnDiskMSnExp"))
    filename <- paste0(path, "/", sapply(fileNames(object), digest::sha1), ".h5")
    if (any(file.exists(filename)))
        stop("File(s) ", paste(filename[file.exists(filename)], collapse = ", "),
             " already exist(s).")
    bpmapply(fileNames(object), filename, FUN = .serialize_msfile_to_hdf5,
             BPPARAM = BPPARAM)
}

readHdf5DiskMSData2 <- function(files, pdata = NULL, msLevel. = NULL,
                               verbose = isMSnbaseVerbose(),
                               centroided. = NA, smoothed. = NA,
                               hdf5path = ".", BPPARAM = bpparam()) {
    obj <- readOnDiskMSData(normalizePath(files), pdata, msLevel.,
                                      verbose, centroided.,
                                      smoothed.)
    hdf5path <- suppressWarnings(normalizePath(hdf5path))
    dir.create(hdf5path, showWarnings = FALSE, recursive = TRUE)
    if (verbose) message("Serialising to hdf5...")
    h5_fls <- serialise_to_hdf52(obj, hdf5path, BPPARAM = BPPARAM)
    obj <- .onDisk2hdf52(obj, h5_fls)
    res <- validHdf5MSnExp2(obj)
    if (is.character(res)) stop(res)
    obj
}

#' Can not simply reuse the `[,OnDiskMSnExp` method as we need to subset also
#' the hdf5file slot accordingly.
#'
#' @author Johannes Rainer
#'
#' @noRd
setMethod("[", "Hdf5MSnExp2", function(x, i, j, ..., drop = TRUE) {
    fnames <- fileNames(x)
    x <- callNextMethod()
    x@hdf5file <- x@hdf5file[match(fileNames(x), fnames)]
    msg <- validHdf5MSnExp2(x)
    if (is.character(msg)) stop(msg)
    else x
})

setMethod("spectrapply", "Hdf5MSnExp2", function(object, FUN = NULL,
                                                 BPPARAM = bpparam(), ...) {
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)
    isOK <- validateFeatureDataForOnDiskMSnExp(fData(object))
    if (!is.null(isOK))
        stop(isOK)
    fDataPerFile <- split.data.frame(fData(object),
                                     f = fData(object)$fileIdx)
    h5_fls <- object@hdf5file
    theQ <- processingQueue(object)
    if (!is.null(FUN))
        theQ <- c(theQ, list(ProcessingStep(FUN, ARGS = list(...))))
    vals <- bpmapply(fDataPerFile, h5_fls,
                     FUN = function(fdata, h5_file, queue) {
                         .apply_processing_queue(
                             .hdf5_read_spectra2(fdata, h5_file),
                             queue)
                     }, MoreArgs = list(queue = theQ),
                     BPPARAM = BPPARAM, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    vals <- unlist(vals, recursive = FALSE)
    vals[rownames(fData(object))]
})

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
#' @param file `character`: the name of the HDF5 file.
#'
#' @param fileNames `character` with the file names of the original data. This
#'     is returned by `fileNames(x)` with `x` being a `Hdf5MSnExp` object.
#'
#' @return list of `Spectrum` objects in the order of the spectra given in
#'     param `fdata`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.hdf5_read_spectra2 <- function(fdata, file) {
    fid <-.Call("_H5Fopen", file, 0L, PACKAGE = "rhdf5")
    on.exit(invisible(.Call("_H5Fclose", fid, PACKAGE = "rhdf5")))
    mzi <- base::lapply(as.character(fdata$spIdx), .h5read_bare, file = fid)
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
