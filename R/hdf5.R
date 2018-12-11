.Hdf5MSnExp <- setClass("Hdf5MSnExp",
                        slots = c(hdf5file = "character",
                                  hdf5handle = "H5IdComponent"),
                        contains = "OnDiskMSnExp",
                        prototype = prototype(
                            new("VersionedBiobase",
                                versions = c(classVersion("MSnExp"),
                                             classVersion("OnDiskMSnExp"),
                                             Hdf5MSnExp = "0.0.1")),
                            spectraProcessingQueue = list(),
                            backend = character()),
                        validity = function(object) {
                            msg <- validMsg(NULL, NULL)
                            if (!length(.Hdf5MSnExp())) {
                                if (!file.exists(object@hdf5file))
                                    msg <- validMsg(msg,
                                                    "hdf5 file is missing.")
                            }
                            if (is.null(msg)) TRUE
                            else msg
                        })


.onDisk2hdf5 <- function(from, filename) {
          to <- .Hdf5MSnExp()
          for (sl in MSnbase:::.slotNames0(from))
              slot(to, sl) <- slot(from, sl)
          to@hdf5file <- filename
          if (validObject(to)) to
      }

serialise_to_hdf5 <- function(object, filename = NULL) {
    stopifnot(inherits(object, "OnDiskMSnExp"))
    if (is.null(filename))
        filename <- paste0(digest::sha1(fileNames(object)), ".h5")
    if (file.exists(filename))
        stop("File ", filename, " already exists.")
    h5 <- rhdf5::H5Fcreate(filename)
    pb <- progress::progress_bar$new(total = length(object))
    for (i in seq_along(fileNames(object))) {
        file_group <- as.character(i)
        file_name  <- fileNames(object)[i]
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
            rhdf5::h5write(.pks, h5, hdfile)
        }
    }
    rhdf5::H5Fclose(h5)
    return(filename)
}

readHdf5DiskMSData <- function(files, pdata = NULL, msLevel. = NULL,
                               verbose = isMSnbaseVerbose(),
                               centroided. = NA, smoothed. = NA,
                               hdf5file = NULL, openHdf5 = TRUE) {
    obj <- MSnbase:::readOnDiskMSData(files, pdata, msLevel.,
                                      verbose, centroided.,
                                      smoothed.)
    if (verbose) message("Serialising to hdf5...")
    hdf5file <- serialise_to_hdf5(obj, hdf5file)
    obj <- .onDisk2hdf5(obj, hdf5file)
    if (openHdf5)
        obj <- hdf5Open(obj)
    if (validObject(obj)) obj
}

hdf5Close <- function(object) {
    if (isHdf5Open(object))
        tryCatch(rhdf5::H5Fclose(object@hdf5handle),
                 error = function(e)
                     stop("Encountered error trying to close ",
                          object@hdf5file,
                          ". Use `h5closeAll()` to force close all hdf5 files."))
    invisible(TRUE)
}


hdf5Open <- function(object) {
    if (!isHdf5Open(object))
        object@hdf5handle <- rhdf5::H5Fopen(object@hdf5file)
    object
}

isHdf5Open <- function(object) {
    stopifnot(inherits(object, "Hdf5MSnExp"))
    rhdf5::H5Iis_valid(object@hdf5handle)
}

setMethod("[[", "Hdf5MSnExp",
          function(x, i, j = "missing", drop = "missing") {
              if (length(i) != 1)
                  stop("subscript out of bounds")
              if (is.character(i))
                  i <- base::match(i, featureNames(x))
              if (any(is.na(i)))
                  stop("subscript out of bounds")
              if (!isHdf5Open(x))
                  x <- hdf5Open(x)
              k <- paste0(fData(x)$fileIdx[[i]], "/",
                          featureNames(x)[[i]])
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
              else MSnbase:::Spectrum2_mz_sorted()
              return(spctr)
          })
