readMSData2 <- function(files,
                        pdata = NULL,
                        msLevel.,
                        verbose = isMSnbaseVerbose(),
                        centroided.,
                        smoothed. = NA) {
    msg <- 'readMSData(..., mode = "onDisk")'
    .Deprecated(paste(strwrap(msg), collapse = "\n"))
    if (missing(msLevel.)) msLevel. <- NULL
    if (missing(centroided.)) centroided. <- NA
    readOnDiskMSData(files = files, pdata = pdata,
                     msLevel. = msLevel., verbose = verbose,
                     centroided. = centroided., smoothed. = smoothed.)
}


readOnDiskMSData <- function(files, pdata, msLevel., verbose,
                             centroided., smoothed.) {
    .testReadMSDataInput(environment())
    stopifnot(is.logical(centroided.))

    ## Creating environment with Spectra objects
    assaydata <- new.env(parent = emptyenv())
    filenams <- filenums <- c()
    fullhd2 <- fullhdorder <- c()
    fullhdordercounter <- 1
    .instrumentInfo <- list()
    ## List eventual limitations
    if (isCdfFile(files)) {
        message("Polarity can not be extracted from netCDF files, please set ",
                "manually the polarity with the 'polarity' method.")
    }
    ## Idea:
    ## o initialize a featureData-data.frame,
    featureDataList <- list()
    ## o for each file, extract header info and put that into featureData
    for (f in files) {
        filen <- match(f, files)
        filenums <- c(filenums, filen)
        filenams <- c(filenams, f)
        ## issue #214: define backend based on file format.
        msdata <- .openMSfile(f)
        .instrumentInfo <- c(.instrumentInfo, list(instrumentInfo(msdata)))
        fullhd <- mzR::header(msdata)
        spidx <- seq_len(nrow(fullhd))
        if (verbose)
            message("Reading ", length(spidx), " spectra from file ",
                    basename(f))
        ## Don't read the individual spectra, just define the names of
        ## the spectra.
        fullhdorder <- c(fullhdorder,
                         formatFileSpectrumNames(fileIds=filen,
                                                 spectrumIds=seq_along(spidx),
                                                 nSpectra=length(spidx),
                                                 nFiles=length(files)))
        ## Extract all Spectrum info from the header and put it into the featureData
        fdData <- fullhd[spidx, , drop = FALSE]
        ## rename totIonCurrent and peaksCount, as detailed in
        ## https://github.com/lgatto/MSnbase/issues/105#issuecomment-229503816
        names(fdData) <- sub("peaksCount", "originalPeaksCount", names(fdData))
        ## Add also:
        ## o fileIdx -> links to fileNames property
        ## o spIdx -> the index of the spectrum in the file.
        fdData <- cbind(fileIdx = rep(filen, nrow(fdData)),
                        spIdx = spidx,
                        smoothed = rep(as.logical(smoothed.), nrow(fdData)),
                        fdData, stringsAsFactors = FALSE)
        if (isCdfFile(f)) {
            ## Add the polarity columns if missing in netCDF
            if (!any(colnames(fdData) == "polarity"))
                fdData <- cbind(fdData, polarity = rep(as.integer(NA),
                                                       nrow(fdData)))
        }
        ## Order the fdData by acquisitionNum to force use of acquisitionNum
        ## as unique ID for the spectrum (issue #103). That way we can use
        ## the spIdx (is the index of the spectrum within the file) for
        ## subsetting and extracting.
        if (!all(sort(fdData$acquisitionNum) == fdData$acquisitionNum))
            warning(paste("Unexpected acquisition number order detected.",
                          "Please contact the maintainers or open an issue",
                          "on https://github.com/lgatto/MSnbase.",
                          sep = "\n")) ## see issue #160
        fdData <- fdData[order(fdData$acquisitionNum), ]
        featureDataList <- c(featureDataList, list(fdData))
        ## Fix for #151; would be nice if we could remove that at some point.
        gc()
        mzR::close(msdata)
        rm(msdata)
    }
    ## new in version 1.9.8
    lockEnvironment(assaydata, bindings = TRUE)
    .cacheEnv <- setCacheEnv(list("assaydata" = assaydata,
                                  "hd" = NULL),
                             level = 0,
                             lock = TRUE)

    ## Create 'MSnProcess' object
    process <- new("MSnProcess",
                   processing = paste0("Data loaded [", date(), "]"),
                   files = files,
                   smoothed = NA)
    ## Create 'fdata' and 'pdata' objects
    if (is.null(pdata)) {
        .pd <- data.frame(sampleNames = basename(files))
        rownames(.pd) <- .pd$sampleNames
        pdata <- new("NAnnotatedDataFrame",
                     data = .pd)
    }
    ## If we've got a featureDataList, use it
    if (length(featureDataList) > 0) {
        fdata <- do.call(rbind, featureDataList)
        fdata <- cbind(fdata, spectrum = 1:nrow(fdata),
                       stringsAsFactors = FALSE)
        ## Setting rownames on the data.frame not on the AnnotatedDataFrame;
        ## did get strange errors otherwise.
        rownames(fdata) <- fullhdorder
        ## Re-order them
        fdata <- fdata[base::sort(fullhdorder), ]
        fdata <- new("AnnotatedDataFrame", data = fdata)
        ## Re-order the features.
        ## fdata <- fdata[ls(assaydata), ]
    } else fdata <- new("AnnotatedDataFrame")

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
    ## Create ProcessingStep if needed.
    ## Create the OnDiskMSnExp object.
    res <- new("OnDiskMSnExp",
               assayData = assaydata,
               phenoData = pdata,
               featureData = fdata,
               processingData = process,
               experimentData = expdata,
               .cache  =  .cacheEnv)
    if (!is.null(msLevel.)) {
        msLevel. <- as.integer(msLevel.)
        res <- filterMsLevel(res, msLevel.)
    }
    if (any(!is.na(centroided.))) {
        if (length(centroided.) == 1) {
            centroided(res) <- centroided.
        } else {
        for (i in seq_along(centroided.))
            centroided(res, msLevel. = i) <- centroided.[i]
        }
    }
    return(res)
}

############################################################
## isCdfFile
##
## Just guessing whether the file is a CDF file based on its ending.
isCdfFile <- function(x) {
    fileEnds <- c("cdf", "nc")
    ## check for endings and and ending followed by a . (e.g. cdf.gz)
    patts <- paste0("\\.", fileEnds, "($|\\.)")
    res <- sapply(patts, function(z) {
        grep(z, x, ignore.case = TRUE)
    })
    return(any(unlist(res)))
}
