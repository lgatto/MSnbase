##' Reads as set of XML-based mass-spectrometry data files and
##' generates an [MSnExp-class] object. This function uses the
##' functionality provided by the `mzR` package to access data and
##' meta data in `mzData`, `mzXML` and `mzML`.
##'
##' When using the `inMemory` mode, the whole MS data is read from
##' file and kept in memory as [Spectrum-class] objects within the
##' [MSnExp-class]'es `assayData` slot.
##'
##' To reduce the memory footpring especially for large MS1 data sets
##' it is also possible to read only selected information from the MS
##' files and fetch the actual spectrum data (i.e. the M/Z and
##' intensity values) only on demand from the original data
##' files. This can be achieved by setting `mode = "onDisk"`. The
##' function returns then an [OnDiskMSnExp-class] object instead of a
##' [MSnExp-class] object.
##'
##' @title Imports mass-spectrometry raw data files as 'MSnExp'
##'     instances.
##' @note `readMSData` uses `normalizePath` to replace relative with
##'     absolute file paths.
##' @aliases readMSData2
##' @md
##' @param files A `character` with file names to be read and parsed.
##' @param pdata An object of class [NAnnotatedDataFrame-class] or
##'     `NULL` (default).
##' @param msLevel. MS level spectra to be read. In `inMemory` mode,
##'     use `1` for MS1 spectra or any larger numeric for MSn
##'     spectra. Default is `2` for `InMemory` mode. `onDisk` mode
##'     supports multiple levels and will, by default, read all the
##'     data.
##' @param verbose Verbosity flag. Default is to use
##'     [isMSnbaseVerbose()].
##' @param centroided. A `logical`, indicating whether spectra are
##'     centroided or not. Default is `NA` in which case the information
##'     is extracted from the raw file (for mzML or mzXML files). In
##'     `onDisk`, it can also be set for different MS levels by a
##'     vector of logicals, where the first element is for MS1, the
##'     second element is for MS2, ... See [OnDiskMSnExp-class] for
##'     an example.
##' @param smoothed. A `logical` indicating whether spectra already
##'     smoothed or not. Default is `NA`.
##' @param cache. Numeric indicating caching level. Default is 0 for
##'     MS1 and 1 MS2 (or higher). Only relevant for `inMemory` mode.
##' @param mode On of `"inMemory"` (default) or `"onDisk"`. The former
##'     loads the raw data in memory, while the latter only generates
##'     the object and the raw data is accessed on disk when
##'     needed. See the *benchmarking* vignette for memory and speed
##'     implications.
##' @return An [MSnExp-class] object for `inMemory` mode and a
##'     [OnDiskMSnExp-class] object for `onDisk` mode.
##' @author Laurent Gatto
##' @seealso [readMgfData()] to read `mgf` peak lists.
##' @examples
##' file <- dir(system.file(package = "MSnbase", dir = "extdata"),
##'             full.name = TRUE,
##'             pattern = "mzXML$")
##' mem <- readMSData(file, mode = "inMemory")
##' mem
##' dsk <- readMSData(file, mode = "onDisk")
##' dsk
readMSData <- function(files, pdata = NULL, msLevel. = NULL,
                       verbose = isMSnbaseVerbose(), centroided. = NA,
                       smoothed. = NA, cache. = 1L,
                       mode = c("inMemory", "onDisk")) {
    mode <- match.arg(mode)
    ## o normalize the file path, i.e. replace relative path with absolute
    ##   path. That fixes possible problems on Windows with SNOW parallel
    ##   processing and also proteowizard problems on unis system with ~ paths.
    files <- normalizePath(files)
    .hasSpecs <- hasSpectra(files)
    suppressWarnings(.hasChroms <- hasChromatograms(files))
    if (any(!.hasSpecs)) {
        msg1 <- paste0("Dropping ", sum(!.hasSpecs),
                       " file(s) without any spectra: ",
                       paste(basename(files[!.hasSpecs]),
                             collapse = ", "), ". ")
        if (all(.hasChroms[!.hasSpecs]))
            msg2 <- "They/it contain(s) chromatograms and can be read with `readSRMData()`."
        else
            msg2 <- paste0("File(s) ",
                           paste(basename(files[!.hasSpecs & .hasChroms]),
                                 collapse = ", "),
                           "contain(s) chromatograms that can be read with `readSRMData`.")
        warning(paste0(msg1, msg2))
    }
    files <- files[.hasSpecs]
    if (!length(files)) {
        process <- new("MSnProcess",
                       processing = paste("No data loaded:", date()))
        if (mode == "inMemory")
            res <- new("MSnExp",
                       processingData = process)
        else res <- new("OnDiskMSnExp",
                        processingData = process)
    } else {
        if (mode == "inMemory") {
            if (is.null(msLevel.)) msLevel. <- 2L
            res <- readInMemMSData(files, pdata = pdata, msLevel. = msLevel.,
                        verbose = verbose, centroided. = centroided.,
                        smoothed. = smoothed., cache. = cache.)
        } else { ## onDisk
            res <- readOnDiskMSData(files = files, pdata = pdata,
                                    msLevel. = msLevel., verbose = verbose,
                                    centroided. = centroided.,
                                    smoothed. = smoothed.)
        }
    }
    res
}


readInMemMSData <- function(files, pdata, msLevel., verbose,
                            centroided., smoothed., cache. = 1) {
    .testReadMSDataInput(environment())
    if (msLevel. == 1) cache. <- 0
    msLevel. <- as.integer(msLevel.)
    ## Creating environment with Spectra objects
    assaydata <- new.env(parent = emptyenv())
    ioncount <- c()
    ioncounter <- 1
    filenams <- filenums <- c()
    fullhd2 <- fullhdorder <- c()
    fullhdordercounter <- 1
    .instrumentInfo <- list()
    ## List eventual limitations
    if (isCdfFile(files)) {
        message("Polarity can not be extracted from netCDF files, please set ",
                "manually the polarity with the 'polarity' method.")
    }
    ## ## Idea:
    ## ## o initialize a featureData-data.frame,
    ## ## o for each file, extract header info and put that into
    ##      featureData;
    for (f in files) {
        filen <- match(f, files)
        filenums <- c(filenums, filen)
        filenams <- c(filenams, f)
        ## issue #214: define backend based on file format.
        msdata <- .openMSfile(f)
        .instrumentInfo <- c(.instrumentInfo, list(instrumentInfo(msdata)))
        fullhd <- mzR::header(msdata)
        ## Issue #325: get centroided information from file, but overwrite if
        ## specified with centroided. parameter.
        if (!is.na(centroided.))
            fullhd$centroided <- as.logical(centroided.)
        spidx <- which(fullhd$msLevel == msLevel.)
        ## increase vectors as needed
        ioncount <- c(ioncount, numeric(length(spidx)))
        fullhdorder <- c(fullhdorder, numeric(length(spidx)))
        if (msLevel. == 1) {
            if (length(spidx) == 0)
                stop("No MS1 spectra in file",f)
            if (verbose) {
                cat("Reading ", length(spidx),
                    " MS1 spectra from file ", basename(f), "\n",
                    sep = "")
                pb <- txtProgressBar(min = 0, max = length(spidx),
                                     style = 3)
            }
            for (i in 1:length(spidx)) {
                if (verbose) setTxtProgressBar(pb, i)
                j <- spidx[i]
                hd <- fullhd[j, ]
                ## Fix missing polarity from netCDF
                pol <- hd$polarity
                if (length(pol) == 0)
                    pol <- NA
                .p <- mzR::peaks(msdata, j)
                sp <- new("Spectrum1",
                          rt = hd$retentionTime,
                          acquisitionNum = as.integer(hd$acquisitionNum),
                          scanIndex = as.integer(hd$seqNum),
                          tic = hd$totIonCurrent,
                          mz = .p[, 1],
                          intensity = .p[, 2],
                          fromFile = as.integer(filen),
                          centroided = as.logical(hd$centroided),
                          smoothed = as.logical(smoothed.),
                          polarity = as.integer(pol))
                ## peaksCount
                ioncount[ioncounter] <- sum(.p[, 2])
                ioncounter <- ioncounter + 1
                .fname <- formatFileSpectrumNames(fileIds=filen,
                                                  spectrumIds=i,
                                                  nSpectra=length(spidx),
                                                  nFiles=length(files))
                assign(.fname, sp, assaydata)
                fullhdorder[fullhdordercounter] <- .fname
                fullhdordercounter <- fullhdordercounter + 1
            }
        } else { ## .msLevel != 1
            if (length(spidx) == 0)
                stop("No MS(n>1) spectra in file", f)
            if (verbose) {
                cat("Reading ", length(spidx), " MS", msLevel.,
                    " spectra from file ", basename(f), "\n",
                    sep = "")
                pb <- txtProgressBar(min = 0, max = length(spidx), style = 3)
            }
            scanNums <- fullhd[fullhd$msLevel == msLevel., "precursorScanNum"]
            if (length(scanNums) != length(spidx))
                stop("Number of spectra and precursor scan number do not match!")
            for (i in 1:length(spidx)) {
                if (verbose) setTxtProgressBar(pb, i)
                j <- spidx[i]
                hd <- fullhd[j, ]
                .p <- mzR::peaks(msdata, j)
                sp <- new("Spectrum2",
                          msLevel = as.integer(hd$msLevel),
                          merged = as.numeric(hd$mergedScan),
                          precScanNum = as.integer(scanNums[i]),
                          precursorMz = hd$precursorMZ,
                          precursorIntensity = hd$precursorIntensity,
                          precursorCharge = as.integer(hd$precursorCharge),
                          collisionEnergy = hd$collisionEnergy,
                          rt = hd$retentionTime,
                          acquisitionNum = as.integer(hd$acquisitionNum),
                          scanIndex = as.integer(hd$seqNum),
                          tic = hd$totIonCurrent,
                          mz = .p[, 1],
                          intensity = .p[, 2],
                          fromFile = as.integer(filen),
                          centroided = as.logical(hd$centroided),
                          smoothed = as.logical(smoothed.),
                          polarity = as.integer(hd$polarity))
                ## peaksCount
                ioncount[ioncounter] <- sum(.p[, 2])
                ioncounter <- ioncounter + 1
                .fname <- formatFileSpectrumNames(fileIds=filen,
                                                  spectrumIds=i,
                                                  nSpectra=length(spidx),
                                                  nFiles=length(files))
                assign(.fname, sp, assaydata)
                fullhdorder[fullhdordercounter] <- .fname
                fullhdordercounter <- fullhdordercounter + 1
            }
        }
        if (cache. >= 1)
            fullhd2 <- rbind(fullhd2, fullhd[spidx, ])
        if (verbose)
            close(pb)
        gc()
        mzR::close(msdata)
        rm(msdata)
    }
    lockEnvironment(assaydata, bindings = TRUE)
    ## cache level 2 yet implemented
    cache. <- testCacheArg(cache., maxCache = 2)
    if (cache. >= 1) {
        fl <- sapply(assaydata, function(x) x@fromFile)
        featnms <- ls(assaydata) ## feature names in final MSnExp
        fl <- fl[featnms] ## reorder file numbers
        stopifnot(all(base::sort(featnms) == base::sort(fullhdorder)))
        fullhdorder <- match(featnms, fullhdorder)
        tmphd <- fullhd2[fullhdorder, ] ## reorder
        ioncount <- ioncount[fullhdorder]
        newhd <- data.frame(fileIdx = fl,
                            retention.time = tmphd$retentionTime,
                            precursor.mz = tmphd$precursorMZ,
                            precursor.intensity = tmphd$precursorIntensity,
                            charge = tmphd$precursorCharge,
                            peaks.count = tmphd$peaksCount,
                            tic = tmphd$totIonCurrent,
                            ionCount = ioncount,
                            ms.level = tmphd$msLevel,
                            acquisition.number = tmphd$acquisitionNum,
                            collision.energy = tmphd$collisionEnergy)
    } else {
        newhd <- NULL ## not used anyway
    }
    .cacheEnv <- setCacheEnv(list("assaydata" = assaydata,
                                  "hd" = newhd),
                             cache., lock = TRUE)
    ## CACHING AS BEEN SUPERSEDED BY THE OnDiskMSnExp IMPLEMENTATION
    ## if cache==2, do not lock assign msdata in .cacheEnv then lock
    ## it and do not close(msdata) above; rm(msdata) is OK

    ## Create 'MSnProcess' object
    process <- new("MSnProcess",
                   processing = paste("Data loaded:", date()),
                   files = files,
                   smoothed = smoothed.)
    ## Create 'fdata' and 'pdata' objects
    nms <- ls(assaydata)
    if (is.null(pdata)) {
        .pd <- data.frame(sampleNames = basename(files))
        rownames(.pd) <- .pd$sampleNames
        pdata <- new("NAnnotatedDataFrame",
                     data = .pd)
    }
    fdata <- new("AnnotatedDataFrame",
                 data = data.frame(
                     spectrum = 1:length(nms),
                     row.names = nms))
    fdata <- fdata[ls(assaydata)] ## reorder features
    ## expriment data slot
    if (length(.instrumentInfo) > 1) {
        cmp <- length(unique(sapply(.instrumentInfo, "[[", 1)))
        if (cmp > 1 & verbose)
            message("According to the instrument information in the files,\n",
                    "the data has been acquired on different instruments!")
        for (nm in names(.instrumentInfo[[1]]))
            .instrumentInfo[[1]][[nm]] <- sapply(.instrumentInfo, "[[", nm)
    }
    expdata <- new("MIAPE",
                   instrumentManufacturer = .instrumentInfo[[1]]$manufacturer,
                   instrumentModel = .instrumentInfo[[1]]$model,
                   ionSource = .instrumentInfo[[1]]$ionisation,
                   analyser = as.character(.instrumentInfo[[1]]$analyzer),
                   detectorType = .instrumentInfo[[1]]$detector)
    ## Create and return 'MSnExp' object
    if (verbose)
        cat("Creating 'MSnExp' object\n")
    toReturn <- new("MSnExp",
                    assayData = assaydata,
                    phenoData = pdata,
                    featureData = fdata,
                    processingData = process,
                    experimentData = expdata,
                    .cache = .cacheEnv)
    return(toReturn)
}

.testReadMSDataInput <- function(e) {
    if (is.numeric(e$msLevel) && !all(e$msLevel > 0))
        stop("msLevel must be an integer > 0.")
    if (length(e$files) < 1)
        stop("At least one MS file is required.")
    if (all(unique(e$files) != e$files))
        stop("Non unique files provided as input. ")
    extensions <- unique(toupper(sub("^.+\\.", "", e$files)))
    if (length(extensions) > 1)
        warning(paste("Reading different file formats in.",
                      "This is untested and you are welcome to try it out.",
                      "Please report back!", sep = "\n"))
    invisible(TRUE)
}
