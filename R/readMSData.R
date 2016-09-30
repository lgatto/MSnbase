readMSData <- function(files,
                       pdata = NULL,
                       msLevel. = 2,
                       verbose = isMSnbaseVerbose(),
                       centroided. = NA,
                       smoothed. = NA,
                       cache. = 1) {
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
    ## ## Idea:
    ## ## o initialize a featureData-data.frame,
    ## ## o for each file, extract header info and put that into
    ##      featureData;
    for (f in files) {
        filen <- match(f, files)
        filenums <- c(filenums, filen)
        filenams <- c(filenams, f)
        msdata <- mzR::openMSfile(f)
        .instrumentInfo <- c(.instrumentInfo, list(instrumentInfo(msdata)))
        fullhd <- mzR::header(msdata)
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
                .p <- mzR::peaks(msdata, j)
                sp <- new("Spectrum1",
                          rt = hd$retentionTime,
                          acquisitionNum = as.integer(hd$acquisitionNum),
                          scanIndex = as.integer(hd$seqNum),
                          tic = hd$totIonCurrent,
                          mz = .p[, 1],
                          intensity = .p[, 2],
                          fromFile = as.integer(filen),
                          centroided = as.logical(centroided.),
                          smoothed = as.logical(smoothed.),
                          polarity = as.integer(hd$polarity))
                ## peaksCount
                ioncount[ioncounter] <- sum(.p[, 2])
                ioncounter <- ioncounter + 1
                .fname <- sprintf(paste0("X%0",
                                         ceiling(log10(length(spidx) + 1L)),
                                         "d.%s"), i, filen)
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
                          centroided = as.logical(centroided.),
                          smoothed = as.logical(smoothed.),
                          polarity = as.integer(hd$polarity))
                ## peaksCount
                ioncount[ioncounter] <- sum(.p[, 2])
                ioncounter <- ioncounter + 1                
                .fname <- sprintf(paste0("X%0",
                                         ceiling(log10(length(spidx) + 1L)),
                                         "d.%s"), i, filen)
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
        newhd <- data.frame(file = fl,
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
                   analyser = .instrumentInfo[[1]]$analyzer,
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
