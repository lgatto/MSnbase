readMSData2 <- function(files,
                        pdata = NULL,
                        msLevel,
                        verbose = TRUE,
                        centroided = FALSE,
                        smoothed = FALSE,
                        removePeaks = 0,
                        clean = FALSE,
                        backend = "disk") {
    .testReadMSDataInput(environment())
    ## Check the backend argument; we're supporting "disk" only for msLevel=1
    ## TODO: add also a trimMz argument.
    ## Creating environment with Spectra objects
    assaydata <- new.env(parent = emptyenv())
    ioncount <- c()
    ioncounter <- 1
    filenams <- filenums <- c()
    fullhd2 <- fullhdorder <- c()
    fullhdordercounter <- 1
    .instrumentInfo <- list()
    ## Idea:
    ## o initialize a featureData-data.frame,
    featureDataList <- list()
    ## o for each file, extract header info and put that into featureData; this might
    ##   be usefull for MS1, but eventually also MS2.
    for (f in files) {
        filen <- match(f, files)
        filenums <- c(filenums, filen)
        filenams <- c(filenams, f)
        msdata <- mzR::openMSfile(f)
        on.exit(close(msdata))
        .instrumentInfo <- c(.instrumentInfo, list(instrumentInfo(msdata)))
        fullhd <- mzR::header(msdata)
        spidx <- seq_len(nrow(fullhd))
        ## increase vectors as needed
        ioncount <- c(ioncount, numeric(length(spidx)))
        ## MS1 level
        if (verbose) {
            message("Reading information of ", length(spidx),
                    " spectra from file ", basename(f), "\n",
                    sep = "")
            pb <- txtProgressBar(min=0, max=length(spidx), style=3)
        }
        ## Don't read the individual spectra, just define the names of the spectra.
        fullhdorder <- c(fullhdorder,
                         sprintf(paste0("X%0",
                                        ceiling(log10(length(spidx) + 1L)),
                                        "d.%s"), 1:length(spidx), filen))
        ## Extract general Spectrum info from the header and put it into the featureData.
        ## This might eventually also be interesting for in-memory MSnExp MS1 data; we might
        ## put this below the if-else.
        ## o acquisitionNum
        ## o polarity
        ## o peaksCount
        ## o totIonCurrent
        ## o retentionTime
        ## o basePeakMZ
        ## o basePeakIntensity
        ## o msLevel
        fdData <- fullhd[spidx, , drop = FALSE]
        ## Add also:
        ## o fileIdx -> links to fileNames property
        ## o spIdx -> the index of the spectrum in the file.
        ## o centroided; the parameter argument.
        fdData <- cbind(fileIdx = filen,
                        spIdx = spidx,
                        centroided = centroided,
                        fdData, stringsAsFactors = FALSE)
        featureDataList <- c(featureDataList, list(fdData))

        if (verbose) setTxtProgressBar(pb, length(spidx))

        ## if (removePeaks > 0)
        ##     sp <- removePeaks(sp, t=removePeaks)
        ## if (clean)
        ##     sp <- clean(sp)

        if (verbose) close(pb)
    }
    ## new in version 1.9.8
    lockEnvironment(assaydata, bindings = TRUE)
    .cacheEnv <- setCacheEnv(list("assaydata" = assaydata,
                                  "hd" = NULL),
                             level = 0,
                             lock = TRUE)
    
    ## and do not close(msdata) above; rm(msdata) is OK
    ## Create 'MSnProcess' object
    process <- new("MSnProcess",
                   processing = paste("Data loaded:",date()),
                   files = files,
                   smoothed = smoothed)
    if (removePeaks > 0) {
        process@processing <- c(process@processing,
                                paste0("Curves <= ", removePeaks, " set to '0': ", date()))
    } else {
        if (clean)
            process@processing <- c(process@processing,
                                    paste("Spectra cleaned: ", date(), sep = ""))
    }
    ## Create 'fdata' and 'pdata' objects
    nms <- ls(assaydata)
    if (is.null(pdata)) {
        .pd <- data.frame(sampleNames = basename(files))
        rownames(.pd) <- .pd$sampleNames
        pdata <- new("NAnnotatedDataFrame",
                     data = .pd)
    }
    ## If we've got the featureDataList, use that one instead; that's for MS1 basically.
    if (length(featureDataList) > 0){
        fdata <- do.call(rbind, featureDataList)
        fdata <- cbind(fdata, spectrum = 1:nrow(fdata), stringsAsFactors = FALSE)
        fdata <- new("AnnotatedDataFrame", data = fdata)
        rownames(fdata) <- fullhdorder
        ## Re-order them
        fdata <- fdata[sort(fullhdorder), ]
        ## Re-order the features.
        ##fdata <- fdata[ls(assaydata), ]
        ## Check if the ordering matches the environment.
        if(!all(ls(assaydata) == rownames(fdata)))
            stop("Ordering of spectra in assayData does not match the order in featureData!")
    }

    ## expriment data slot
    if (length(.instrumentInfo) > 1) {
        cmp <- sapply(.instrumentInfo[-1], function(x) identical(x, .instrumentInfo[[1]]))
        if (!all(cmp)) {
            warning("According to the instrument information in the files, the data has been acquired on different instruments!")
            .instrumentInfo[[1]] <-
                list(manufacturer = paste(sapply(.instrumentInfo, "[[", "manufacturer"), collapse = ", "),
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
    ## Create ProcessingStep if needed.
    queue <- list()
    if (removePeaks > 0)
        queue <- c(queue,
                   list(ProcessingStep(FUN = "removePeaks",
                                       ARGS = list(t = removePeaks))))
    if (clean)
        queue <- c(queue,
                   list(ProcessingStep(FUN = "clean")))
    ## Create the OnDiskMSnExp object.
    if (verbose)
        cat("Creating 'MSnExp' object\n")
    res <- new("OnDiskMSnExp",
               assayData = assaydata,
               phenoData = pdata,
               featureData = fdata,
               processingData = process,
               experimentData = expdata,
               spectraProcessingQueue = queue,
               .cache  =  .cacheEnv)
    if (!missing(msLevel)) {
        msLevel <- as.integer(msLevel)
        res <- res[fData(res)$msLevel %in% msLevel]
    }
    return(res)
}

