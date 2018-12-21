setClassUnion("H5IdComponentOrNULL", c("H5IdComponent", "NULL"))

##' @title The `Hdf5MSnExp` Class for MS Data And Meta-Data
##' @description The `Hdf5MSnExp` class encapsulates data and
##'     meta-data for mass spectrometry experiments like the `MSnExp`
##'     and `OnDiskMSnExp` classes. `Hdf5MSnExp` implements, like
##'     `OnDiskMSnExp`, an *on disk* model, where the raw data (M/Z
##'     and intensities) are stored on disk (rather than in memory,
##'     like `MSnExp` data) using, as the name implies, the hdf5 file
##'     format. See [`MSnExp`] and [`OnDiskMSnExp`] for a more general
##'     description of the data, slots and operations.
##'
##'     Object of the class are currently created with the
##'     `readHdf5DiskMSData` function. Later, this backend will be
##'     included in the main `readMSData` function.
##'
##' @slot hdf5file `character(1)` containing the hdf5 filename.
##' @slot hdf5handle A `H5IdComponent` or `NULL`.
##'
##' @rdname Hdf5MSnExp-class
##' @aliases Hdf5MSnExp
##' @seealso See [`MSnExp`] and [`OnDiskMSnExp`] classes.
##' @author Laurent Gatto
##' @md
##' @examples
##' f <- msdata::proteomics(pattern = "MS3TMT11", full.names = TRUE)
##' x <- readHdf5DiskMSData(f)
##' x
##'
##' ## automatically generated hdf5 file
##' hdf5FileName(x)
##' hdf5FileName(x) <- "myhdf5data.h5"
##'
##' x[[1]]
##' x[[2]]
##' x[[10]]
##' filterMsLevel(x, 3L)[[1]]
##'
##' ## clean up session
##' file.remove(hdf5FileName(x))
##' validHdf5MSnExp(x) ## not valid anymore
##' rm(x)
.Hdf5MSnExp <- setClass("Hdf5MSnExp",
                        slots = c(hdf5file = "character",
                                  hdf5handle = "H5IdComponentOrNULL"),
                        contains = "OnDiskMSnExp",
                        prototype = prototype(
                            new("VersionedBiobase",
                                versions = c(classVersion("MSnExp"),
                                             classVersion("OnDiskMSnExp"),
                                             Hdf5MSnExp = "0.0.1")),
                            spectraProcessingQueue = list(),
                            backend = character()))

validHdf5MSnExp <- function(object) {
    if (!length(.Hdf5MSnExp())) {
        if (!file.exists(object@hdf5file))
            stop("hdf5 file is missing.")
    }
    TRUE
}

.onDisk2hdf5 <- function(from, filename) {
          to <- .Hdf5MSnExp()
          for (sl in MSnbase:::.slotNames0(from))
              slot(to, sl) <- slot(from, sl)
          to@hdf5file <- filename
          to
      }

serialise_to_hdf5 <- function(object, filename = NULL) {
    stopifnot(inherits(object, "OnDiskMSnExp"))
    if (is.null(filename))
        filename <- paste0(digest::sha1(fileNames(object)), ".h5")
    if (file.exists(filename))
        stop("File ", filename, " already exists.")
    h5 <- rhdf5::H5Fcreate(filename)
    pb <- progress::progress_bar$new(total = length(object))
    comp_level <- .hdf5_compression_level()
    for (i in seq_along(fileNames(object))) {
        file_name  <- fileNames(object)[i]
        file_group <- .hdf5_group_name(file_name)
        stopifnot(rhdf5::h5createGroup(h5, file_group))
        fns <- featureNames(filterFile(object, i))
        fh <- openMSfile(file_name)
        pks <- peaks(fh)
        close(fh)
        for (j in seq_along(fns)) {
            pb$tick()
            fn <- fns[j]
            hdfile <- paste0(file_group, "/", fn)
            .pks <- pks[[j]]
            colnames(.pks) <- c("mz", "intensity")
            rhdf5::h5write(.pks, h5, hdfile, level = comp_level)
        }
    }
    rhdf5::H5Fclose(h5)
    return(filename)
}

##' @rdname Hdf5MSnExp-class
readHdf5DiskMSData <- function(files, pdata = NULL, msLevel. = NULL,
                               verbose = isMSnbaseVerbose(),
                               centroided. = NA, smoothed. = NA,
                               hdf5file = NULL, openHdf5 = TRUE) {
    obj <- MSnbase:::readOnDiskMSData(normalizePath(files), pdata, msLevel.,
                                      verbose, centroided.,
                                      smoothed.)
    if (verbose) message("Serialising to hdf5...")
    hdf5file <- serialise_to_hdf5(obj, hdf5file)
    obj <- .onDisk2hdf5(obj, hdf5file)
    if (openHdf5)
        obj <- hdf5Open(obj)
    if (validHdf5MSnExp(obj)) obj
}

##' The `hdf5Close` and `hdf5Open` function respectively close and
##' open the connection to the `Hdf5MSnExp` object hdf5 file. The both
##' return the input object with the updated hdf5 file
##' handle. `isHdf5Open` test if the hdf5 file is open and retuns a
##' `logical(1)`.
##'
##' @param object An instance of class `Hdf5MSnExp`.
##' @rdname Hdf5MSnExp-class
hdf5Close <- function(object) {
    if (isHdf5Open(object))
        tryCatch(rhdf5::H5Fclose(object@hdf5handle),
                 error = function(e)
                     stop("Encountered error trying to close ",
                          object@hdf5file,
                          ". Use `h5closeAll()` to force close all hdf5 files."))
    invisible(TRUE)
}

##' @rdname Hdf5MSnExp-class
hdf5Open <- function(object) {
    if (!isHdf5Open(object))
        object@hdf5handle <- rhdf5::H5Fopen(object@hdf5file)
    object
}

##' @rdname Hdf5MSnExp-class
isHdf5Open <- function(object) {
    stopifnot(inherits(object, "Hdf5MSnExp"))
    if (is.null(object@hdf5handle))
        return(FALSE)
    rhdf5::H5Iis_valid(object@hdf5handle)
}

##' @rdname Hdf5MSnExp-class
hdf5FileName <- function(object) {
    stopifnot(inherits(object, "Hdf5MSnExp"))
    object@hdf5file
}

##' @rdname Hdf5MSnExp-class
`hdf5FileName<-` <- function(object, value) {
    stopifnot(inherits(object, "Hdf5MSnExp"))
    stopifnot(inherits(value, "character"))
    if (length(value) != 1)
        stop("Please provide a single file name")
    if (isHdf5Open(object))
        object <- hdf5Close(object)
    stopifnot(file.rename(hdf5FileName(object), value))
    object@hdf5file <- value
    if (validHdf5MSnExp(object))
        return(object)
}



setMethod("[[", "Hdf5MSnExp",
          function(x, i, j = "missing", drop = "missing") {
              if (!isHdf5Open(x))
                  x <- hdf5Open(x)
              if (length(i) != 1)
                  stop("subscript out of bounds")
              if (is.character(i))
                  i <- base::match(i, featureNames(x))
              if (any(is.na(i)))
                  stop("subscript out of bounds")
              if (!isHdf5Open(x))
                  x <- hdf5Open(x)
              k <- paste0(.hdf5_group_name(fileNames(x)[fData(x)$fileIdx[[i]]]),
                                           "/", featureNames(x)[[i]])
              rw <- rhdf5::h5read(x@hdf5handle, k)
              if (msLevel(x)[i] == 1L)
                  spctr <- MSnbase:::Spectrum1_mz_sorted(
                                         rt = rtime(x)[[i]],
                                         acquisitionNum = acquisitionNum(x)[[i]],
                                         scanIndex = scanIndex(x)[[i]],
                                         tic = tic(x)[[i]],
                                         mz = rw[, 1],
                                         intensity = rw[, 2],
                                         fromFile = fromFile(x)[[i]],
                                         centroided = centroided(x)[[i]],
                                         smoothed = smoothed(x)[[i]],
                                         polarity = polarity(x)[[i]])
              else
                  spctr <- MSnbase:::Spectrum2_mz_sorted(
                                         msLevel = msLevel(x)[[i]],
                                         rt = rtime(x)[[i]],
                                         acquisitionNum = acquisitionNum(x)[[i]],
                                         scanIndex = scanIndex(x)[[i]],
                                         tic = tic(x)[[i]],
                                         mz = rw[, 1],
                                         intensity = rw[, 2],
                                         fromFile = fromFile(x)[[i]],
                                         centroided = centroided(x)[[i]],
                                         smoothed = smoothed(x)[[i]],
                                         polarity = polarity(x)[[i]],
                                         precScanNum = precScanNum(x)[[i]],
                                         precursorMz = precursorIntensity(x)[[i]],
                                         precursorIntensity = precursorIntensity(x)[[i]],
                                         precursorCharge = precursorCharge(x)[[i]],
                                         collisionEnergy = collisionEnergy(x)[[i]])
              if (validObject(spctr)) return(spctr)
          })

.hdf5_group_name <- function(x) {
    vapply(x, digest::sha1, character(1), USE.NAMES = FALSE)
}

#' This is Mike Smith's function (issue #395) that directly access the data
#' without validation and file checking.
#'
#' @param handle `H5IdComponent` as stored in the `@hdf5handle` slot of the
#'     `Hdf5MSnExp`.
#'
#' @param name `character` defining the data set to be read.
#'
#' @return the imported data set (in most cases a `matrix`).
#'
#' @author Mike Smith
#'
#' @md
#'
#' @noRd
.h5read_raw <- function(file, name = "") {
    did <- .Call("_H5Dopen", file@ID, name, NULL, PACKAGE = "rhdf5")
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
#' @param fileh `H5IdComponent` as stored in the `@hdf5handle` slot of the
#'     `Hdf5MSnExp`.
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
.read_spectra_hdf5 <- function(fdata, fileh, fileNames) {
    group_names <- MSnbase:::.hdf5_group_name(fileNames[fdata$fileIdx])
    k <- paste0(group_names, "/", rownames(fdata))
    mzi <- lapply(k, .h5read_optimised, file = fileh)
    res <- vector("list", nrow(fdata))
    names(res) <- rownames(fdata)
    ms1 <- which(fdata$msLevel == 1)
    n_peaks <- base::lengths(mzi, use.names = FALSE) / 2
    if (length(ms1)) {
        mzi_ms1 <- do.call(rbind, mzi[ms1])
        res[ms1] <- MSnbase:::Spectra1_mz_sorted(
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
        res[msn] <- MSnbase:::Spectra2_mz_sorted(
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
