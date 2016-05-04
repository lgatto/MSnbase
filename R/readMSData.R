.testReadMSDataInput <- function(e) {
    if (!all(e$msLevel > 0))
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

## readMSData2 <- function(files,
##                        pdata = NULL,
##                        msLevel = 2,
##                        verbose = TRUE,
##                        centroided = FALSE,
##                        smoothed = FALSE,
##                        removePeaks = 0,
##                        clean = FALSE,
##                        cache = 2) {
## ## - deal with multiple files
##   if (removePeaks > 0)
##     stop("Can't remove peaks when cache level is > 1")
##   if (clean > 0)
##     stop("Can't clean peaks when cache level is > 1")
##   MSnbase:::.testReadMSDataInput(environment())
##   if (msLevel == 1) ## cache currently only works for MS2 level data
##     cache <- 0
##   msLevel <- as.integer(msLevel)
##   mzr <- lapply(files, mzR::openMSfile)
##   .instrumentInfo <- lapply(mzr, instrumentInfo)
##   names(.instrumentInfo) <- names(mzr) <- seq_along(mzr)
##   fullhd <- do.call(rbind, lapply(mzr, header))
##   ioncount <- unlist(lapply(mzr, function(xx)
##                             sapply(peaks(xx), function(pk) sum(pk[, 2]))))
##   .n <- sapply(mzr, length)
##   featnms <- lapply(.n, function(n) paste0("X", 1:n))
##   featnms <- unlist(lapply(seq_along(featnms),
##                            function(x) {
##                              paste(featnms[[x]], x, sep=".")
##                            }))
##   rownames(fullhd) <- featnms
##   fromfile <- unlist(mapply(rep, 1:length(.n), .n, SIMPLIFY = FALSE)) ## to be saved to cache
##   newhd <- data.frame(file = fromfile,
##                       retention.time = fullhd$retentionTime,
##                       precursor.mz = fullhd$precursorMZ,
##                       precursor.intensity = fullhd$precursorIntensity,
##                       charge = fullhd$precursorCharge,
##                       peaks.count = fullhd$peaksCount,
##                       tic = fullhd$totIonCurrent,
##                       ionCount = ioncount,
##                       ms.level = fullhd$msLevel,
##                       acquisition.number = fullhd$acquisitionNum,
##                       collision.energy = fullhd$collisionEnergy)
##   ## if (verbose)
##   ##   message("Caching...")
##   ## .cacheEnv <- MSnbase:::setCacheEnv(list("assaydata" = assaydata,
##   ##                                         "hd" = newhd),
##   ##                                    cache, lock = TRUE)
##   ## Create 'MSnProcess' object
##   process <- new("MSnProcess",
##                  processing = c(
##                    paste("Data loaded:", date()),
##                    paste("Cache level ", cache)),
##                  files = files,
##                  smoothed = smoothed)
##   if (is.null(pdata)) {
##     .pd <- data.frame(sampleNames = files,
##                       fileNumbers = seq_along(files))
##     pdata <- new("NAnnotatedDataFrame",
##                  data = .pd)
##   }
##   fdata <- new("AnnotatedDataFrame",
##                data=data.frame(
##                  spectrum=1:nrow(fullhd) ,
##                  row.names=featnms))
##   ## expriment data slot
##   if (length(.instrumentInfo) > 1) {
##     cmp <- sapply(.instrumentInfo[-1], function(x) identical(x, .instrumentInfo[[1]]))
##     if (!all(cmp)) {
##       warning("According to the instrument information in the files, the data has been acquired on different instruments!")
##       .instrumentInfo[[1]] <- list(manufacturer = paste(sapply(.instrumentInfo, "[[", "manufacturer"), collapse = ", "),
##                                    model = paste(sapply(.instrumentInfo, "[[", "model"), collapse = ", "),
##                                    ionisation = paste(sapply(.instrumentInfo, "[[", "ionisation"), collapse = ", "),
##                                    analyzer = paste(sapply(.instrumentInfo, "[[", "analyzer"), collapse = ", "),
##                                    detector = paste(sapply(.instrumentInfo, "[[", "detector"), collapse = ", "))
##     }
##   }
##   expdata <- new("MIAPE",
##                  instrumentManufacturer = .instrumentInfo[[1]]$manufacturer,
##                  instrumentModel = .instrumentInfo[[1]]$model,
##                  ionSource = .instrumentInfo[[1]]$ionisation,
##                  analyser = .instrumentInfo[[1]]$analyzer,
##                  detectorType = .instrumentInfo[[1]]$detector)
##   ## Create and return 'MSnExp' object
##   if (verbose)
##     cat("Creating 'MSnExp' object\n")
##   toReturn <- new("MSnExp",
##                   assayData = list2env(mzr),
##                   phenoData = pdata,
##                   featureData = fdata,
##                   processingData = process,
##                   experimentData = expdata)
##                   ## .cache = .cacheEnv)
##   if (validObject(toReturn))
##     return(toReturn)
## }

readMSData <- function(files,
                       pdata = NULL,
                       msLevel = 2,
                       verbose = TRUE,
                       centroided = FALSE,
                       smoothed = FALSE,
                       removePeaks = 0,
                       clean = FALSE,
                       cache = 1) {
    .testReadMSDataInput(environment())
    ## TODO: add also a trimMz argument.
    if (msLevel == 1) ## cache currently only works for MS2 level data
        cache <- 0
    msLevel <- as.integer(msLevel)
    ## Creating environment with Spectra objects
    assaydata <- new.env(parent=emptyenv())
    ioncount <- c()
    ioncounter <- 1
    filenams <- filenums <- c()
    fullhd2 <- fullhdorder <- c()
    fullhdordercounter <- 1
    .instrumentInfo <- list()
    for (f in files) {
        filen <- match(f, files)
        filenums <- c(filenums, filen)
        filenams <- c(filenams, f)
        msdata <- mzR::openMSfile(f)
        .instrumentInfo <- c(.instrumentInfo, list(instrumentInfo(msdata)))
        fullhd <- mzR::header(msdata)
        ifelse(msLevel == 1, ## later, > 1 level
               spidx <- which(fullhd$msLevel == 1),
               spidx <- which(fullhd$msLevel > 1))
        ## increase vectors as needed
        ioncount <- c(ioncount, numeric(length(spidx)))
        fullhdorder <- c(fullhdorder, numeric(length(spidx)))
        ## MS1 level
        if (msLevel == 1) {
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
                sp <- new("Spectrum1",
                          peaksCount = hd$peaksCount,
                          rt = hd$retentionTime,
                          acquisitionNum = hd$acquisitionNum,
                          mz = mzR::peaks(msdata,j)[,1],
                          intensity = mzR::peaks(msdata,j)[,2],
                          fromFile = filen,
                          tic = hd$totIonCurrent,
                          centroided = centroided)
                ioncount[ioncounter] <- sum(mzR::peaks(msdata,j)[,2])
                ioncounter <- ioncounter + 1
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
            if (length(spidx)==0)
                stop("No MS(n>1) spectra in file",f)
            if (verbose) {
                cat("Reading ",length(spidx)," MS2 spectra from file ",
                    basename(f),"\n",sep="")
                pb <- txtProgressBar(min=0,max=length(spidx),style=3)
            }
            ## this was fullhd$acquisitionNum -- check/wrong
            ## ms1scanNums <- MSnbase:::getBins(fullhd$msLevel[spidx])
            scanNums <- fullhd[fullhd$msLevel == 2, "precursorScanNum"]
            if (length(scanNums)!=length(spidx))
                stop("Number of spectra and precursor scan number do not match!")
            for (i in 1:length(spidx)) {
                if (verbose) setTxtProgressBar(pb, i)
                j <- spidx[i]
                hd <- fullhd[j,]
                .p <- mzR::peaks(msdata,j)
                sp <- new("Spectrum2",
                          precScanNum = as.integer(scanNums[i]),
                          precursorMz = hd$precursorMZ,
                          precursorIntensity = hd$precursorIntensity,
                          precursorCharge = hd$precursorCharge,
                          collisionEnergy = hd$collisionEnergy,
                          tic = hd$totIonCurrent,
                          peaksCount = hd$peaksCount,
                          rt = hd$retentionTime,
                          acquisitionNum = hd$acquisitionNum,
                          mz = .p[,1],
                          intensity = .p[,2],
                          fromFile = filen,
                          centroided = centroided)
                ioncount[ioncounter] <- sum(.p[,2])
                ioncounter <- ioncounter + 1
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
        }
        if (cache >= 1)
            fullhd2 <- rbind(fullhd2, fullhd[spidx, ])
        if (verbose)
            close(pb)
        gc() ## could this help with Error in function (x): no function to return from, jumping to top level)...
        mzR::close(msdata) ## DO NOT CLOSE IF CACHE LEVEL >= 2
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
    ## if cache==2, do not lock,
    ## assign msdata in .cacheEnv
    ## then lock it
    ## and do not close(msdata) above; rm(msdata) is OK
    ## Create 'MSnProcess' object
    process <- new("MSnProcess",
                   processing = paste("Data loaded:",date()),
                   files = files,
                   smoothed = smoothed)
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
        .pd <- data.frame(sampleNames = 1)
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
        if (!all(cmp)) {
            warning("According to the instrument information in the files, the data has been acquired on different instruments!")
            .instrumentInfo[[1]] <- list(manufacturer = paste(sapply(.instrumentInfo, "[[", "manufacturer"), collapse = ", "),
                                         model = paste(sapply(.instrumentInfo, "[[", "model"), collapse = ", "),
                                         ionisation = paste(sapply(.instrumentInfo, "[[", "ionisation"), collapse = ", "),
                                         analyzer = paste(sapply(.instrumentInfo, "[[", "analyzer"), collapse = ", "),
                                         detector = paste(sapply(.instrumentInfo, "[[", "detector"), collapse = ", "))
        }
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
    if (validObject(toReturn))
        return(toReturn)
}

