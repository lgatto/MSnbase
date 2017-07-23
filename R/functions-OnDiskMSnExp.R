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

#' @title Apply a function to spectra loaded from a single file
#' 
#' @description This function creates `Spectrum1` and `Spectrum2` objects for
#'     the specicied spectra in one file and applies the provided function to
#'     each of them.
#'
#' @note Using the C constructor that takes all values at once and creates a
#'     list of `Spectrum1` objects, applies processing steps, applies the
#'     provided function and returns its results - or the list of
#'     `Spectrum1` objects if `APPLYFUN = NULL`.
#'
#'     Note: enforces ordering of M/Z-intensity pairs by M/Z.
#' 
#' @param fData: either a full `data.frame` (returned by `fData(OnDiskMSnExp))`
#'     or a sub-set forspecific spectra. The data.frame should ONLY
#'     CONTAIN VALUES FOR SPECTRA OF ONE FILE!
#'
#' @param filenames: `fileNames(object)` with `object` being an `OnDiskMSnExp`.
#'
#' @param queue: `object@spectraProcessingQueue`; if `lenght > 0` all
#'     processing steps will be applied to the created spectrum
#'     objects.
#' 
#' @param APPLYFUN: the function to be applied to the spectrum objects
#'     (such as `ionCount` etc). If `NULL` the function returns the list of
#'     spectrum objects.
#' 
#' @param fastLoad: `logical(1)` whether reading the spectras' header data
#'     should be omitted prior to retrieving the data (i.e. skip the
#'     `mzR::header` call before calling `mzR::peaks`. The former call might be
#'     required on some systems (so far macOS) for some files.
#' 
#' @param ...: additional arguments for the APPLYFUN
#'
#' @return `list` with either spectrum objects or the results of the function
#'     provided with argument `APPLYFUN`.
#' 
#' @noRd
#'
#' @md
.applyFun2SpectraOfFileMulti <- function(fData, filenames,
                                         queue = NULL,
                                         APPLYFUN = NULL,
                                         fastLoad = TRUE,
                                         ...) {
    suppressPackageStartupMessages(
        require(MSnbase, quietly = TRUE)
    )
    verbose. <- isMSnbaseVerbose()
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are required!")
    filename <- filenames[fData[1, "fileIdx"]]
    ## issue #214: define backend based on file format.
    fileh <- .openMSfile(filename)
    msLevel1 <- which(fData$msLevel == 1)
    msLevelN <- which(fData$msLevel > 1)
    ## Reading the header for the selecte spectra. This is to avoid getting
    ## "memory not mapped" errors when reading mz and intensity values from
    ## certain mzML files (issue #170). Since this problem seems to be absent
    ## on linux and Windows systems we allow the user to disable it.
    ## Also we are just reading the header for the last spectrum since that
    ## seems to fix it too.
    if (!fastLoad)
        hd_spectra <- mzR::header(fileh, max(fData$spIdx))
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
        rm(res2)
    }
    mzR::close(fileh)
    rm(fileh)
    ## Intermediate #151 fix. Performance-wise would be nice to get rid of this.
    ## gc()
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
    ## Ensure that ordering is the same than in fData:
    res[match(rownames(fData), names(res))]
}

#' @description Less memory demanding version of `applyFun2SpectraOfFileMulti`.
#'
#' @details This function loops through the rows of the `data.frame` passed
#'     with argument `fData`, loads the mz and intensity values of the
#'     respective spectrum from the original file, creates a `Spectrum` object,
#'     executes eventual *lazy evaluation* steps (argument `queue`) and
#'     executes, if provided, the `APPLYFUN` to the `Spectrum`.
#'
#' @note The function uses the *C* constructor for `Spectrum` objects that
#'     enforces also ordering of M/Z-intensity pairs by M/Z.
#' 
#' @param fData: either a full `data.frame` (returned by `fData(OnDiskMSnExp))`
#'     or a sub-set forspecific spectra. The data.frame should ONLY
#'     CONTAIN VALUES FOR SPECTRA OF ONE FILE!
#'
#' @param filenames: `fileNames(object)` with `object` being an `OnDiskMSnExp`.
#'
#' @param queue: `object@spectraProcessingQueue`; if `lenght > 0` all
#'     processing steps will be applied to the created spectrum
#'     objects.
#' 
#' @param APPLYFUN: the function to be applied to the spectrum objects
#'     (such as `ionCount` etc). If `NULL` the function returns the list of
#'     spectrum objects.
#' 
#' @param fastLoad: `logical(1)` whether reading the spectras' header data
#'     should be omitted prior to retrieving the data (i.e. skip the
#'     `mzR::header` call before calling `mzR::peaks`. The former call might be
#'     required on some systems (so far macOS) for some files.
#' 
#' @param ...: additional arguments for the APPLYFUN
#'
#' @return `list` with either spectrum objects or the results of the function
#'     provided with argument `APPLYFUN`.
#'
#' @author Johannes Rainer
#' 
#' @noRd
#' @md
.applyFun2IndividualSpectraOfFile <- function(fData, filenames,
                                                  queue = NULL,
                                                  APPLYFUN = NULL,
                                                  fastLoad = TRUE,
                                                  ...) {
    suppressPackageStartupMessages(
        require(MSnbase, quietly = TRUE)
    )
    verbose. <- isMSnbaseVerbose()
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are required!")
    if (length(queue) > 0) {
        if (verbose.) {
            message("Apply lazy processing step(s):")
            for (j in 1:length(queue))
                message(" o '", queue[[j]]@FUN, "' with ",
                        length(queue[[j]]@ARGS), " argument(s).")
        }
    }
    filename <- filenames[fData[1, "fileIdx"]]
    ## issue #214: define backend based on file format.
    fileh <- .openMSfile(filename)
    ## Reading the header for the selecte spectra. This is to avoid getting
    ## "memory not mapped" errors when reading mz and intensity values from
    ## certain mzML files (issue #170). Since this problem seems to be absent
    ## on linux and Windows systems we allow the user to disable it.
    ## Also we are just reading the header for the last spectrum since that
    ## seems to fix it too.
    if (!fastLoad)
        hd_spectra <- mzR::header(fileh, max(fData$spIdx))
    n_rows <- nrow(fData)
    res <- vector("list", n_rows)
    for (i in 1:n_rows) {
        pks <- mzR::peaks(fileh, fData$spIdx[i])
        if (fData$msLevel[i] == 1) {
            sp <- Spectrum1_mz_sorted(rt = fData$retentionTime[i],
                                      acquisitionNum = fData$acquisitionNum[i],
                                      scanIndex = fData$spIdx[i],
                                      tic = fData$totIonCurrent[i],
                                      mz = pks[, 1],
                                      intensity = pks[, 2],
                                      fromFile = fData$fileIdx[i],
                                      centroided = fData$centroided[i],
                                      smoothed = fData$smoothed[i],
                                      polarity = fData$polarity[i])
        } else {
            sp <- Spectrum2_mz_sorted(msLevel = fData$msLevel[i],
                                      rt = fData$retentionTime[i],
                                      acquisitionNum = fData$acquisitionNum[i],
                                      scanIndex = fData$spIdx[i],
                                      tic = fData$totIonCurrent[i],
                                      mz = pks[, 1],
                                      intensity = pks[, 2],
                                      fromFile = fData$fileIdx[i],
                                      centroided = fData$centroided[i],
                                      smoothed = fData$smoothed[i],
                                      polarity = fData$polarity[i],
                                      merged = fData$mergedScan[i],
                                      precScanNum = fData$precursorScanNum[i],
                                      precursorMz = fData$precursorMZ[i],
                                      precursorIntensity = fData$precursorIntensity[i],
                                      precursorCharge = fData$precursorCharge[i],
                                      collisionEnergy = fData$collisionEnergy[i])
        }
        ## And now go through the processing queue - if not empty...
        if (length(queue) > 0) {
            for (pStep in queue) {
                sp <- executeProcessingStep(pStep, sp)
            }
        }
        ## Apply the function, if provided
        if (!is.null(APPLYFUN))
            sp <- APPLYFUN(sp, ...)
        res[[i]] <- sp
    }
    names(res) <- rownames(fData)
    mzR::close(fileh)
    rm(fileh)
    ## Intermediate #151 fix. Performance-wise would be nice to get rid of this.
    ## gc()
    res
}


## Extract chromatogram(s).
#' @description This function extracts chromatograms efficiently for multiple
#'     rt and mz ranges by loading the data per file only once and performing
#'     the mz subsetting on the already loaded Spectrum1 classes.
#'
#' @note x has to be an MSnExp or OnDiskMSnExp object. The function is optimized
#'     for OnDiskMSnExp objects such that extraction is performed in parallel
#'     for each file.
#' 
#' @param rt \code{matrix} with two columns and number of rows corresponding to
#'     the number of ranges to extract. If the number of columns of the matrix
#'     is not equal to 2, \code{range} is called on each row.
#'
#' @param mz \code{matrix} with two columns and number of rows corresponding to
#'     the number of ranges to extract. nrow of rt and mz have to match. If the
#'     number of columns of the matrix is not equal to 2, \code{range} is
#'     called on each row.
#'
#' @param x OnDiskMSnExp object from which to extract the chromatograms.
#'
#' @param missingValue value to be used as intensity if no signal was measured
#'     for a given rt.
#'
#' @param msLevel \code{integer(1)} ensuring that the chromatogram is extracted
#'     only for a specified MS level.
#' 
#' @return A \code{matrix} with the \code{Chromatogram} objects with rows
#'     corresponding to ranges and columns to files/samples. \code{result[, 1]}
#'     will thus return a \code{list} of \code{Chromatogram} objects for the
#'     first sample/file, while \code{result[1, ]} returns a \code{list} of
#'     \code{Chromatogram} objects for the same rt/mz range for all files.
#'
#' @author Johannes Rainer
#'
#' @noRd
.extractMultipleChromatograms <- function(x, rt, mz, aggregationFun = "sum",
                                          BPPARAM = bpparam(),
                                          missingValue = NA_real_,
                                          msLevel = 1L) {
    missingValue <- as.numeric(missingValue)
    if (!any(.SUPPORTED_AGG_FUN_CHROM == aggregationFun))
        stop("'aggregationFun' should be one of ",
             paste0("'", .SUPPORTED_AGG_FUN_CHROM, "'", collapse = ", "))
    ## Ensure we're working on one MS level only!
    fns <- fileNames(x)
    x <- filterMsLevel(x, msLevel)
    if (length(x) == 0) {
        warning("No MS ", msLevel, " data present.")
        return(matrix(ncol = length(fns), nrow = 0))
    }
    if (missing(rt))
        rt <- matrix(c(-Inf, Inf), nrow = 1)
    if (missing(mz))
        mz <- matrix(c(-Inf, Inf), nrow = 1)
    ## Calculate the range for each row in rt 
    if (ncol(rt) != 2)
        rt <- t(apply(rt, MARGIN = 1, range))
    ## Replicate if nrow rt is 1 to match nrow of mz.
    if (nrow(rt) == 1)
        rt <- matrix(rep(rt, nrow(mz)), ncol = 2, byrow = TRUE)
    if (ncol(mz) != 2)
        mz <- t(apply(mz, MARGIN = 1, range))
    if (nrow(mz) == 1)
        mz <- matrix(rep(mz, nrow(rt)), ncol = 2, byrow = TRUE)
    if (nrow(rt) != nrow(mz))
        stop("dimensions of 'rt' and 'mz' have to match")

    ## Identify indices of all spectra that are within the rt ranges.
    rtimes <- rtime(x)
    keep_idx <- unlist(apply(rt, MARGIN = 1, function(z)
        which(rtimes >= z[1] & rtimes <= z[2])), use.names = FALSE)
    keep_idx <- sort(unique(as.integer(keep_idx)))
    if (length(keep_idx) == 0)
        return(matrix(ncol = length(fns), nrow = 0))
    ## 1) Subset x keeping all spectra that fall into any of the provided rt
    ##    ranges.
    subs <- x[keep_idx]

    ## 2) Call the final subsetting on each file separately.
    subs_by_file <- splitByFile(subs, f = factor(seq_along(fileNames(subs))))
    suppressWarnings(
        res <- bpmapply(
            subs_by_file,
            match(fileNames(subs), fns),
            FUN = function(cur_sample, cur_file, rtm, mzm, aggFun) {
                ## Load all spectra for that file. applies also any proc steps
                sps <- spectra(cur_sample)
                ## Related to issue #229: can we avoid getting all spectra and
                ## just return the intensity values for each spectrum instead?
                rts <- rtime(cur_sample)
                cur_res <- vector("list", nrow(rtm))
                ## Loop through rt and mz.
                for (i in 1:nrow(rtm)) {
                    ## - Select all spectra within that range and call a
                    ##   function on them that does first filterMz and then
                    ##   aggregate the values per spectrum.
                    in_rt <- rts >= rtm[i, 1] & rts <= rtm[i, 2]
                    ## Return an empty Chromatogram if there is no spectrum/scan
                    ## within the retention time range.
                    if (!any(in_rt)) {
                        cur_res[[i]] <- Chromatogram(
                            filterMz = mzm[i, ],
                            fromFile = as.integer(cur_file),
                            aggregationFun = aggFun)
                        next
                    }
                    cur_sps <- lapply(
                        sps[in_rt],
                        function(spct, filter_mz, aggFun) {
                            spct <- filterMz(spct, filter_mz)
                            ## Now aggregate the values.
                            if (!spct@peaksCount)
                                return(c(NA_real_, NA_real_, missingValue))
                            c(range(spct@mz, na.rm = TRUE, finite = TRUE),
                              do.call(aggFun, list(spct@intensity,
                                                   na.rm = TRUE)))
                        }, filter_mz = mzm[i, ], aggFun = aggFun)
                    ## Now build the Chromatogram class.
                    allVals <- unlist(cur_sps, use.names = FALSE)
                    idx <- seq(3, length(allVals), by = 3)
                    ## Or should we drop the names completely?
                    ints <- allVals[idx]
                    names(ints) <- names(cur_sps)
                    ## If no measurement is non-NA, still report the NAs and
                    ## use the filter mz as mz. We hence return a
                    ## Chromatogram with retention times of all spectra
                    ## within the mz/rt slice, but with all intensity
                    ## values being NA
                    mz_range <- mzm[i, ]
                    if (!all(is.na(ints)))
                        mz_range <- range(allVals[-idx], na.rm = TRUE,
                                          finite = TRUE)
                    cur_res[[i]] <- Chromatogram(
                        rtime = rts[in_rt],
                        intensity = ints,
                        mz = mz_range,
                        filterMz = mzm[i, ],
                        fromFile = as.integer(cur_file),
                        aggregationFun = aggFun)
                }
                cur_res
            }, MoreArgs = list(rtm = rt, mzm = mz, aggFun = aggregationFun),
            BPPARAM = BPPARAM, SIMPLIFY = FALSE)
    )
    ## Ensure that the lists have the same length than there are samples!
    fromF <- base::match(fileNames(subs), fns)

    ## Jumping here if we have a file with no spectra within the rt range.
    ## If we've got some files in which we don't have any signal in any range,
    ## fill it with empty Chromatograms. This ensures that the result has
    ## ALWAYS the same length/number of columns than there are samples.
    if (length(res) != length(fns)) {
        res_all_files <- vector(mode = "list", length = length(fns))
        res_all_files[fromF] <- res
        empties <- which(lengths(res_all_files) == 0)
        ## fill these
        for (i in empties) {
            empty_list <- vector(mode = "list", length = nrow(rt))
            for(j in 1:nrow(rt)) {
                empty_list[[j]] <- Chromatogram(filterMz = mz[j, ],
                                                fromFile = as.integer(i),
                                                aggregationFun = aggregationFun)
            }
            res_all_files[[i]] <- empty_list
        }
        res <- res_all_files
    }
    do.call(cbind, res)
}
