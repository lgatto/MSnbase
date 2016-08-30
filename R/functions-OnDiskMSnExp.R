############################################################
## subsetFeatureDataBy
##
## Convenience function to subset a OnDiskMSnExp featureData based on
## provided subsetting criteria:
## o index: subset by numeric index, logical or character. Character
##    indices are "forwarded" to argument "name". index refers to the
##    rows in the feature data.
## o scanIdx: a numeric is expected, specifying the scan index. If a
##    character vector is provided, it is "forwarded" to
##    argument "name".
## o scanIdxCol: the column of the featureData containing the scan indices.
## o name: a character vector, matching is performed using the row names
##    of the fd data.frame.
## o rtlim: a numeric of length 2 specifying the retention time window
##    from which spectra should be extracted.
## o msLevel: a numeric of MS level.
## Open questions:
## 1) Should we re-calculate the spIdx or the seqNum if we're subsetting?
## Note and comments:
## o At present filtering by multiple criteria (e.g. index and name) is not
##    possible.
subsetFeatureDataBy <- function(fd, index = NULL, scanIdx = NULL,
                                scanIdxCol = "acquisitionNum", name = NULL,
                                rtlim = NULL, msLevel = NULL) {
    valMsg <- validateFeatureDataForOnDiskMSnExp(fd)
    if (is.character(valMsg))
        stop(valMsg)
    ## o index
    if (length(index) > 0) {
        if (is.logical(index)) {
            if (length(index) != nrow(fd))
                stop("If 'index' is a logical vector its length has to match",
                     " the number of",
                     " rows of the featureData!")
            index <- which(index)
        }
        if (is.numeric(index)) {
            gotIt <- index %in% 1:nrow(fd)
            if (!any(gotIt))
                stop("Provided indices are outside of the allowed range.")
            if (any(!gotIt))
                warning("Some of the provided indices are outside of the",
                        " allowed range.")
            index <- index[gotIt]
            return(fd[index, , drop = FALSE])
        }
        if (is.character(index))
            name <- index
    }
    ## o scanIdx
    if (length(scanIdx) > 0) {
        if(is.numeric(scanIdx)){
            gotIt <- scanIdx %in% fd[, scanIdxCol]
            if(!any(gotIt))
                stop("None of the provided scan indices are available!")
            if(!all(gotIt))
                warning("Some of the provided scan indices are not available.")
            return(fd[which(fd[, scanIdxCol] %in% scanIdx), , drop=FALSE])
        }
        if (is.character(scanIdx))
            name <- scanIdx
    }
    ## o name: subset by name, match to rownames.
    if (length(name) > 0) {
        gotIt <- name %in% rownames(fd)
        if (!any(gotIt))
            stop("None of the provided names found.")
        if (!all(gotIt))
            warning("Some of the provided names do not match featureData rownames.")
        name <- name[gotIt]
        return(fd[name, , drop=FALSE])
    }
    ## o rtlim: subset by retention time range.
    if (length(rtlim > 0)){
        if (length(rtlim) > 2 | !is.numeric(rtlim))
            stop("Argument 'rtlim' has to be a numeric vector of length 2 specifying",
                 " the retention time window (range).")
        gotIt <- which(fd$retentionTime >= rtlim[1] & fd$retentionTime <= rtlim[2])
        if (length(gotIt) == 0)
            stop("No spectrum within the specified retention time window.")
        fd <- fd[gotIt, , drop=FALSE]
    }
    ## o msLevel
    fd <- fd[fd$msLevel %in% msLevel, , drop = FALSE]
    if (nrow(fd) == 0)
        warning("No spectra with the specified MS level present.")
    return(fd)
}

## Returns either NULL or a character string.
validateFeatureDataForOnDiskMSnExp <- function(x) {
    ## Testing if we've got all the required columns! See issue 105
    ## for a discussion about originalTotIonCurrent and
    ## originalPeaksCount.
    reqCols <- c("fileIdx", "spIdx", "acquisitionNum",
                 "retentionTime", "polarity", "msLevel",
                 "totIonCurrent", "originalPeaksCount",
                 "centroided")
    NotPresent <- reqCols[!(reqCols %in% colnames(x))]
    if (length(NotPresent) > 0)
        return(paste0("Required columns: ",
                      paste(sQuote(NotPresent), collapse = ","),
                      " not present in featureData!"))
    return(NULL)
}

############################################################
## validateOnDiskMSnExp
##
## The validation method that might be called manually. In addition to the
## validate function called by validObject this ensures also that all
## spectra objects are valid and thus re-reads the raw data.
## o mzTolerance: allowed tolerance in comparison of spectras' M/Z ranges
##    with featureData's lowMZ and highMZ.
validateOnDiskMSnExp <- function(object, mzTolerance=1e-6) {
    ## First call the basic validity.
    validObject(object)
    ## Now check validity of the spectra; if one non-valid object is found we stop.
    tmp <- spectrapply(object, FUN = function(z) {
        res <- validObject(z)
        if (is(res, "character"))
            stop(res)
    })
    ## Spectra M/Z range: spectra M/Z ranges should be within lowMZ and highMZ.
    ## Get the mzrange by spectrum.
    res <- spectrapply(object, FUN = function(x) {
        return(range(mz(x), na.rm=TRUE))
    })
    ## Loop through the spectra and check that mzrange is within
    ## lowMZ and highMZ
    fd <- fData(object)
    emptyMz <- 0
    for (i in seq_len(length(object))) {
        ## Check if lowMZ or highMZ are 0; in that case skip and throw
        ## a warning, see https://github.com/lgatto/MSnbase/issues/153
        if (fd$lowMZ[i] == 0L & fd$highMZ[i] == 0L) {
            emptyMz <- emptyMz + 1
        } else {
            hgeq <- isTRUE(base::all.equal(res[[i]][2], fd$highMZ[i],
                                           tolerance = mzTolerance))
            lweq <- isTRUE(base::all.equal(res[[i]][1], fd$lowMZ[i],
                                           tolerance = mzTolerance))
            if (!(hgeq & lweq))
                stop("M/Z values of spectrum ", i, " don't match ",
                     "lowMZ and highMZ.")
        }
    }
    if (emptyMz > 0)
        warning("Header value lowMZ and highMZ is 0 for ",
                emptyMz, " spectra.")
    
    ## Check that msLevel of spectra matches msLevel of featureData.
    if (any(fData(object)$msLevel != sapply(spectra(object), msLevel)))
        stop("msLevel in featureData does not match msLevel from the spectra!")
    return(TRUE)
}

############################################################
## spectrapply
##
## That's the main method to apply functions to the object's spectra, or
## to just return a list with the spectra, if FUN is empty.
## Parallel processing by file can be enabled using BPPARAM.
spectrapply <- function(object, FUN = NULL,
                        BPPARAM = bpparam(), ...) {
    if (!is(object, "OnDiskMSnExp"))
        stop("'object' is expected to be an 'OnDiskMSnExp' object!")
    ## Check if we would do better with serial processing:
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)
    isOK <- validateFeatureDataForOnDiskMSnExp(fData(object))
    if (!is.null(isOK))
        stop(isOK)
    fDataPerFile <- base::split(fData(object),
                                f = fData(object)$fileIdx)
    fNames <- fileNames(object)
    theQ <- processingQueue(object)
    vals <- bplapply(fDataPerFile,
                     FUN = .applyFun2SpectraOfFileMulti,
                     filenames = fNames,
                     queue = theQ,
                     APPLYFUN = FUN,
                     BPPARAM = BPPARAM,
                     ...)
    names(vals) <- NULL
    vals <- unlist(vals, recursive = FALSE)
    return(vals[rownames(fData(object))])
}

############################################################
## precursorValue_OnDiskMSnExp
##
## Returns requested information from the featureData, ensuring
## that data for MS1 is set to NA. Throws an error if the
## object contains only MS1 data.
precursorValue_OnDiskMSnExp <- function(object, column) {
    ## Throw an error if we've got only MS1:
    if (all(unique(msLevel(object)) == 1))
        stop("This experiment contains MS1 spectra.")
    ps <- fData(object)[, column]
    names(ps) <- featureNames(object)
    ## Replacing values for MS 1 with NA.
    ps[msLevel(object) == 1] <- NA
    return(ps)
}

## Using the C constructor that takes all values at once and creates a
## list of Spectrum1 objects, applies processing steps, applies the
## provided function and returns its results - or the list of
## Spectrum1 objects if APPLYFUN = NULL.
## Note: enforces ordering of M/Z-intensity pairs by M/Z.
## Arguments:
## o fData: either a full data.frame (returned by fData(OnDiskMSnExp))
##   or a sub-set forspecific spectra. The data.frame should ONLY
##   CONTAIN VALUES FOR SPECTRA OF ONE FILE!
## o filenames: fileNames(object)
## o queue: object@spectraProcessingQueue; if lenght > 0 all
##   processing steps will be applied to the created Spectrum1
##   objects.
## o APPLYFUN: the function to be applied to the Spectrum1 objects
##   (such as ionCount etc).  If NULL the function returns the list of
##   Spectrum1 objects.
## o ...: additional arguments for the APPLYFUN
.applyFun2SpectraOfFileMulti <- function(fData, filenames,
                                         queue = NULL,
                                         APPLYFUN = NULL,
                                         ...) {
    suppressPackageStartupMessages(
        require(MSnbase, quietly = TRUE)
    )
    verbose. <- isMSnbaseVerbose()
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are required!")
    filename <- filenames[fData[1, "fileIdx"]]
    ## Open the file.
    fileh <- mzR::openMSfile(filename)
    ## hd <- header(fileh)
    on.exit(expr = mzR::close(fileh))
    ## Intermediate #151 fix. Performance-wise would be nice to get rid of this.
    on.exit(expr = gc(), add = TRUE)
    msLevel1 <- which(fData$msLevel == 1)
    msLevelN <- which(fData$msLevel > 1)
    ## Process MS1 and MSn separately
    if (length(msLevel1) >= 1) {
        ms1fd <- fData[msLevel1, , drop = FALSE]
        ## Reading all of the data in "one go". According to issue
        ## #103 we should use acquisitionNum, not spectrum idx.
        ## See issue #118 for an explanation of the match
        ## allSpect <- mzR::peaks(fileh,
        ##                        match(ms1fd$acquisitionNum, hd$acquisitionNum))
        ## spIdx is the index of the spectra in the provided header, thus we
        ## don't require to load the header and do the matching.
        allSpect <- mzR::peaks(fileh, ms1fd$spIdx)
        ## If we have more than one spectrum the peaks function returns a list.
        if (is(allSpect, "list")) {
            nValues <- base::lengths(allSpect, use.names = FALSE) / 2
            allSpect <- do.call(rbind, allSpect)
        } else {
            ## otherwise it's a matrix, e.g. if only a single scan
            ## index was provided.
            nValues <- nrow(allSpect)
        }
        ## Call the C-constructor to create a list of Spectrum1
        ## objects.
        ## Sorting of M/Z values as discussed in issue #135
        ## Benchmarks for this: issue #136
        res <- Spectra1_mz_sorted(peaksCount = nValues,
                                  rt = ms1fd$retentionTime,
                                  acquisitionNum = ms1fd$acquisitionNum,
                                  scanIndex = ms1fd$spIdx,
                                  tic = ms1fd$totIonCurrent,
                                  mz = allSpect[, 1],
                                  intensity = allSpect[, 2],
                                  fromFile = ms1fd$fileIdx,
                                  centroided = ms1fd$centroided,
                                  smoothed = ms1fd$smoothed,
                                  polarity = ms1fd$polarity,
                                  nvalues = nValues)
        names(res) <- rownames(ms1fd)
    } else {
        res <- list()
    }
    if (length(msLevelN) >= 1) {
        msnfd <- fData[msLevelN, , drop = FALSE]
        ## Reading all of the data in "one go".
        ## See issue #118 for an explanation of the match/no match.
        ## allSpect <- mzR::peaks(fileh,
        ##                        match(msnfd$acquisitionNum, hd$acquisitionNum))
        allSpect <- mzR::peaks(fileh, msnfd$spIdx)
        ## If we have more than one spectrum the peaks function returns a list.
        if (is(allSpect, "list")) {
            nValues <- as.integer(base::lengths(allSpect, use.names = FALSE) / 2)
            allSpect <- do.call(rbind, allSpect)
        } else {
            ## otherwise it's a matrix, e.g. if only a single scan
            ## index was provided.
            nValues <- nrow(allSpect)
        }
        ## Call the C-constructor to create a list of Spectrum2
        ## objects.
        ## Sorting of M/Z values as discussed in issue #135
        res2 <- Spectra2_mz_sorted(msLevel = msnfd$msLevel,
                                   peaksCount = nValues,
                                   rt = msnfd$retentionTime,
                                   acquisitionNum = msnfd$acquisitionNum,
                                   scanIndex = msnfd$spIdx,
                                   tic = msnfd$totIonCurrent,
                                   mz = allSpect[, 1],
                                   intensity = allSpect[, 2],
                                   fromFile = msnfd$fileIdx,
                                   centroided = msnfd$centroided,
                                   smoothed = msnfd$smoothed,
                                   polarity = msnfd$polarity,
                                   merged = msnfd$mergedScan,
                                   precScanNum = msnfd$precursorScanNum,
                                   precursorMz = msnfd$precursorMZ,
                                   precursorIntensity = msnfd$precursorIntensity,
                                   precursorCharge = msnfd$precursorCharge,
                                   collisionEnergy = msnfd$collisionEnergy,
                                   nvalues = nValues)
        names(res2) <- rownames(msnfd)
        res <- c(res, res2)
    }
    ## Ensure that ordering is the same than in fData:
    res <- res[match(rownames(fData), names(res))]
    ## If we have a non-empty queue, we might want to execute that too.
    if (!is.null(APPLYFUN) | length(queue) > 0){
        if (length(queue) > 0) {
            if (verbose.) {
                message("Apply lazy processing step(s):")
                for (j in 1:length(queue))
                    message(" o '", queue[[j]]@FUN, "' with ",
                            length(queue[[j]]@ARGS), " argument(s).")
            }
        }
        res <- lapply(res, FUN = function(z, theQ, APPLF, ...){
            ## Apply the processing steps.
            if (length(theQ) > 0) {
                for (pStep in theQ) {
                    z <- executeProcessingStep(pStep, z)
                }
            }
            if (is.null(APPLF)) {
                return(z)
            } else {
                return(do.call(APPLF, args = c(list(z), ...)))
            }
        }, theQ = queue, APPLF = APPLYFUN, ...)
    }
    return(res)
}


## Same as above, but using a for loop and the C-constructor for individual
## Spectrum objects.
.applyFun2SpectraOfFileSingle <- function(fData, filenames,
                                         queue = NULL,
                                         APPLYFUN = NULL,
                                         ...) {
    suppressPackageStartupMessages(
        require(MSnbase, quietly = TRUE)
    )
    verbose. <- isMSnbaseVerbose()
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are required!")
    filename <- filenames[fData[1, "fileIdx"]]
    ## Open the file.
    fileh <- mzR::openMSfile(filename)
    hd <- header(fileh)
    on.exit(expr = mzR::close(fileh))
    msLevel1 <- which(fData$msLevel == 1)
    msLevelN <- which(fData$msLevel > 1)
    ## Process MS1 and MSn separately
    if (length(msLevel1) >= 1) {
        ms1fd <- fData[msLevel1, , drop = FALSE]
        ## Reading all of the data in "one go". According to issue
        ## #103 we should use acquisitionNum, not spectrum idx.
        ## See issue #118 for an explanation of the match
        allSpect <- mzR::peaks(fileh,
                               match(ms1fd$acquisitionNum, hd$acquisitionNum))
        ## If we have more than one spectrum the peaks function returns a list.
        if (!is(allSpect, "list"))
            allSpect <- list(allSpect)
        ## Do it with a for loop.
        res <- vector("list", nrow(ms1fd))
        for (i in 1:nrow(ms1fd)) {
            currentMat <- allSpect[[i]]
            o <- order(currentMat[, 1], method = "radix")
            currentMat <- currentMat[o, ]
            res[[i]] <- Spectrum1(peaksCount = nrow(currentMat),
                                  scanIndex = ms1fd[i, "spIdx"],
                                  rt = ms1fd[i, "retentionTime"],
                                  acquisitionNum = ms1fd[i, "acquisitionNum"],
                                  mz = currentMat[, 1],
                                  intensity = currentMat[, 2],
                                  centroided = ms1fd[i, "centroided"],
                                  smoothed = ms1fd[i, "smoothed"],
                                  fromFile = ms1fd[i, "fileIdx"],
                                  polarity = ms1fd[i, "polarity"],
                                  tic = ms1fd[i, "totIonCurrent"]
                                  )
        }
        names(res) <- rownames(ms1fd)
    } else {
        res <- list()
    }
    if (length(msLevelN) >= 1) {
        msnfd <- fData[msLevelN, , drop = FALSE]
        ## Reading all of the data in "one go".
        ## See issue #118 for an explanation of the match
        allSpect <- mzR::peaks(fileh,
                               match(msnfd$acquisitionNum, hd$acquisitionNum))
        ## If we have more than one spectrum the peaks function returns a list.
        if (!is(allSpect, "list"))
            allSpect <- list(allSpect)
        ## Do it with a for loop.
        res2 <- vector("list", nrow(msnfd))
        for (i in 1:nrow(msnfd)) {
            currentMat <- allSpect[[i]]
            o <- order(currentMat[, 1], method = "radix")
            currentMat <- currentMat[o, ]
            res2[[i]] <- Spectrum2(peaksCount = nrow(currentMat),
                                   scanIndex = msnfd[i, "spIdx"],
                                   rt = msnfd[i, "retentionTime"],
                                   acquisitionNum = msnfd[i, "acquisitionNum"],
                                   mz = currentMat[, 1],
                                   intensity = currentMat[, 2],
                                   centroided = msnfd[i, "centroided"],
                                   smoothed = msnfd[i, "smoothed"],
                                   fromFile = msnfd[i, "fileIdx"],
                                   polarity = msnfd[i, "polarity"],
                                   tic = msnfd[i, "totIonCurrent"],
                                   msLevel = msnfd[i, "msLevel"],
                                   merged = msnfd[i, "mergedScan"],
                                   precScanNum = msnfd[i, "precursorScanNum"],
                                   precursorMz = msnfd[i, "precursorMZ"],
                                   precursorIntensity = msnfd[i, "precursorIntensity"],
                                   precursorCharge = msnfd[i, "precursorCharge"],
                                   collisionEnergy = msnfd[i, "collisionEnergy"]
                                   )
        }
        names(res2) <- rownames(msnfd)
        res <- c(res, res2)
    }
    ## Ensure that ordering is the same than in fData:
    res <- res[match(rownames(fData), names(res))]
    ## If we have a non-empty queue, we might want to execute that too.
    if (!is.null(APPLYFUN) | length(queue) > 0){
        if (length(queue) > 0) {
            if (verbose.) {
                message("Apply lazy processing step(s):")
                for (j in 1:length(queue))
                    message(" o '", queue[[j]]@FUN, "' with ",
                            length(queue[[j]]@ARGS), " argument(s).")
            }
        }
        res <- lapply(res, FUN = function(z, theQ, APPLF, ...){
            ## Apply the processing steps.
            if (length(theQ) > 0) {
                for (pStep in theQ) {
                    z <- executeProcessingStep(pStep, z)
                }
            }
            if (is.null(APPLF)) {
                return(z)
            } else {
                return(do.call(APPLF, args = c(list(z), ...)))
            }
        }, theQ = queue, APPLF = APPLYFUN, ...)
    }
    return(res)
}

