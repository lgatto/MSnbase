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
                      "Please report back!", sep="\n"))
    invisible(TRUE)
}

readMSData <- function(files,
                       pdata = NULL,
                       msLevel. = 2,
                       verbose = TRUE,
                       centroided. = FALSE,
                       smoothed. = FALSE,
                       removePeaks = 0,
                       clean = FALSE,
                       cache = 1) {
    .testReadMSDataInput(environment())
    ## TODO: add also a trimMz argument.
    if (msLevel. == 1) ## cache currently only works for MS2 level data
        cache <- 0
    msLevel. <- as.integer(msLevel.)
    ## Creating environment with Spectra objects
    assaydata <- new.env(parent = emptyenv())
    ioncount <- c()
    filenams <- filenums <- c()
    fullhd2 <- fullhdorder <- c()
    fullhdordercounter <- 1
    .instrumentInfo <- list()
    ## ## Idea:
    ## ## o initialize a featureData-data.frame,
    ## featureDataList <- list()
    ## ## o for each file, extract header info and put that into featureData; this might
    ## ##   be usefull for MS1, but eventually also MS2.
    for (f in files) {
        filen <- match(f, files)
        filenums <- c(filenums, filen)
        filenams <- c(filenams, f)
        msdata <- mzR::openMSfile(f)
        on.exit(mzR::close(msdata)) 
        .instrumentInfo <- c(.instrumentInfo, list(instrumentInfo(msdata)))
        fullhd <- mzR::header(msdata)
        ifelse(msLevel. == 1, ## later, > 1 level
               spidx <- which(fullhd$msLevel == 1),
               spidx <- which(fullhd$msLevel > 1))
        ## increase vectors as needed
        ioncount <- c(ioncount, numeric(length(spidx)))
        fullhdorder <- c(fullhdorder, numeric(length(spidx)))
        ## MS1 level
        if (msLevel. == 1) {
            if (length(spidx) == 0)
                stop("No MS1 spectra in file",f)
            if (verbose) {
                cat("Reading ", length(spidx), " MS1 spectra from file ",
                    basename(f),"\n",sep="")
                pb <- txtProgressBar(min=0, max=length(spidx), style=3)
            }
            for (i in 1:length(spidx)) {
                if (verbose) setTxtProgressBar(pb, i)
                j <- spidx[i]
                hd <- fullhd[j,]
                pks <- mzR::peaks(msdata, j)
                sp <- new("Spectrum1",
                          rt = hd$retentionTime,
                          acquisitionNum = hd$acquisitionNum,
                          mz = pks[, 1],
                          scanIndex = hd$seqNum,
                          tic = hd$totIonCurrent,
                          intensity = pks[, 2],
                          fromFile = filen,
                          polarity = hd$polarity,
                          centroided = centroided.)
                if (removePeaks > 0)
                    sp <- removePeaks(sp, t=removePeaks)
                if (clean)
                    sp <- clean(sp)
                .fname <- sprintf(paste0("X%0",
                                         ceiling(log10(length(spidx) + 1L)),
                                         "d.%s"), i, filen)
                assign(.fname, sp, assaydata)
                fullhdorder[fullhdordercounter] <- .fname
                fullhdordercounter <- fullhdordercounter + 1
            }
        } else { ## MS>2 levels
            if (length(spidx) == 0)
                stop("No MS(n>1) spectra in file",f)
            if (verbose) {
                cat("Reading ", length(spidx),
                    " MS2 spectra from file ", basename(f), "\n",
                    sep = "")
                pb <- txtProgressBar(min = 0, max = length(spidx), style = 3)
            }
            scanNums <- fullhd[fullhd$msLevel == 2, "precursorScanNum"]
            if (length(scanNums) != length(spidx))
                stop("Number of spectra and precursor scan number do not match!")
            for (i in 1:length(spidx)) {
                if (verbose) setTxtProgressBar(pb, i)
                j <- spidx[i]
                hd <- fullhd[j,]
                .p <- mzR::peaks(msdata, j)
                sp <- new("Spectrum2",
                          scanIndex = as.integer(hd$seqNum),
                          merged = as.numeric(hd$mergedScan),
                          polarity = hd$polarity,
                          precScanNum = as.integer(scanNums[i]),
                          precursorMz = hd$precursorMZ,
                          precursorIntensity = hd$precursorIntensity,
                          precursorCharge = hd$precursorCharge,
                          collisionEnergy = hd$collisionEnergy,
                          rt = hd$retentionTime,
                          acquisitionNum = hd$acquisitionNum,
                          mz = .p[, 1],
                          tic = hd$totIonCurrent,
                          intensity = .p[, 2],
                          fromFile = filen,
                          centroided = centroided.)
                if (removePeaks > 0)
                    sp <- removePeaks(sp, t = removePeaks)
                if (clean)
                    sp <- clean(sp)
                .fname <- sprintf(paste0("X%0",
                                         ceiling(log10(length(spidx) + 1L)),
                                         "d.%s"), i, filen)
                assign(.fname, sp, assaydata)
                fullhdorder[fullhdordercounter] <- .fname
                fullhdordercounter <- fullhdordercounter + 1
            }
        }
        if (cache >= 1)
            fullhd2 <- rbind(fullhd2, fullhd[spidx, ])
        if (verbose)
            close(pb)
        rm(msdata)
    }
    ## new in version 1.9.8
    lockEnvironment(assaydata, bindings = TRUE)
    ## cache level 2 yet implemented
    cache <- testCacheArg(cache, maxCache = 2)
    if (cache >= 1) {
        ## results sometimes in:
        ##  Error in function (x)  : attempt to apply non-function
        ## fl <- sapply(assaydata, fromFile)
        fl <- sapply(assaydata, function(x) x@fromFile)
        featnms <- ls(assaydata) ## feature names in final MSnExp
        fl <- fl[featnms] ## reorder file numbers
        stopifnot(all(sort(featnms) == sort(fullhdorder)))
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
    if (verbose)
        message("Caching...")
    .cacheEnv <- setCacheEnv(list("assaydata" = assaydata,
                                  "hd" = newhd),
                             cache, lock = TRUE)
    ## CACHING AS BEEN SUPERSEDED BY THE OnDiskMSnExp IMPLEMENTATION
    ## if cache==2, do not lock assign msdata in .cacheEnv then lock
    ## it and do not close(msdata) above; rm(msdata) is OK
    
    ## Create 'MSnProcess' object
    process <- new("MSnProcess",
                   processing = paste("Data loaded:",date()),
                   files = files,
                   smoothed = smoothed.)
    if (removePeaks > 0) {
        process@processing <- c(process@processing,
                                paste("Curves <= ", removePeaks, " set to '0': ", date(), sep=""))
    } else {
        if (clean)
            process@processing <- c(process@processing,
                                    paste("Spectra cleaned: ", date(), sep=""))
    }
    ## Create 'fdata' and 'pdata' objects
    nms <- ls(assaydata)
    if (is.null(pdata)) {
        .pd <- data.frame(sampleNames = basename(files))
        rownames(.pd) <- .pd$sampleNames
        pdata <- new("NAnnotatedDataFrame",
                     data = .pd)
    }
    fdata <- new("AnnotatedDataFrame",
                 data=data.frame(
                     spectrum=1:length(nms),
                     row.names=nms))
    fdata <- fdata[ls(assaydata)] ## reorder features
    ## expriment data slot
    if (length(.instrumentInfo) > 1) {
        cmp <- sapply(.instrumentInfo[-1], function(x) identical(x, .instrumentInfo[[1]]))
#### START ORIGINAL CODE
        ## if (!all(cmp)) {
        ##     warning("According to the instrument information in the files, the data has been acquired on different instruments!")
        ##     .instrumentInfo[[1]] <- list(manufacturer = paste(sapply(.instrumentInfo, "[[", "manufacturer"), collapse = ", "),
        ##                                  model = paste(sapply(.instrumentInfo, "[[", "model"), collapse = ", "),
        ##                                  ionisation = paste(sapply(.instrumentInfo, "[[", "ionisation"), collapse = ", "),
        ##                                  analyzer = paste(sapply(.instrumentInfo, "[[", "analyzer"), collapse = ", "),
        ##                                  detector = paste(sapply(.instrumentInfo, "[[", "detector"), collapse = ", "))
        ## }
#### END ORIGINAL CODE
#### START FROM readMSData2
        for (nm in names(.instrumentInfo[[1]]))
            .instrumentInfo[[1]][[nm]] <- sapply(.instrumentInfo, "[[", nm)

#### END FROM readMSData2
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
    ## if (validObject(toReturn))  ## validity checks are already performed with "new", no need to perform twice.
    return(toReturn)
}

