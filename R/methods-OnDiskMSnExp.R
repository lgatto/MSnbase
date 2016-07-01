############################################################
## Methods for OnDiskMSnExp objects.

############################################################
## initialize
##
setMethod("initialize",
          signature(.Object="OnDiskMSnExp"),
          function(.Object, ...){
              ## if(!any(names(list(...)) == "onDisk")){
              ##     .Object@onDisk=TRUE
              ## }
              callNextMethod()
          })

setMethod("header",
          c("OnDiskMSnExp", "missing"),
          function(object) {
              .Deprecated("fData")
              fData(object)
          })

setMethod("header",
          c("OnDiskMSnExp", "numeric"),
          function(object, scans) {
              .Deprecated("fData")
              fData(object)[scans, ]
          })


############################################################
## processingQueue
setMethod("processingQueue", "OnDiskMSnExp", function(object){
    return(object@spectraProcessingQueue)
})

############################################################
## msLevel
##
## Extract the msLevel info for all spectra in an OnDiskMSnExp
## object. In contrast to the MSnExp we're not getting that
## from the individual Spectrum objects in @assayData, but
## from the featureData.
setMethod("msLevel", "OnDiskMSnExp", function(object){
    msl <- featureData(object)$msLevel
    ## That should not happen!
    if(is.null(msl))
        stop("The 'OnDiskMSnExp' does not have msLevel information in the featureData!")
    names(msl) <- featureNames(object)
    return(msl)
})

############################################################
## fromFile
##
## Extract the index from which file the spectra in the data set
## derive from.
setMethod("fromFile", "OnDiskMSnExp", function(object){
    fidx <- fData(object)$fileIdx
    names(fidx) <- featureNames(object)
    return(fidx)
})

############################################################
## length
##
## Get the length, i.e. the number of spectra we've got (from the
## featureData).
setMethod("length", "OnDiskMSnExp", function(x){
    return(nrow(fData(x)))
})

############################################################
## scanIndex
##
## Get the scan index for each spectrum in each file. We're extracting
## that from the featureData.
setMethod("scanIndex", "OnDiskMSnExp", function(object){
    scIdx <- fData(object)$spIdx
    names(scIdx) <- featureNames(object)
    return(scIdx)
})

############################################################
## acquisitionNum
##
## Get the acquisition number for each spectrum in each file. We're extracting
## that from the featureData.
setMethod("acquisitionNum", "OnDiskMSnExp", function(object){
    aIdx <- fData(object)$acquisitionNum
    names(aIdx) <- featureNames(object)
    return(aIdx)
})

############################################################
## centroided
##
## Getter/setter for the centroided information; extracting this from
## the featureData.
setMethod("centroided","OnDiskMSnExp",
          function(object){
              val <- fData(object)$centroided
              names(val) <- featureNames(object)
              return(val)
          })
setReplaceMethod("centroided", signature(object="OnDiskMSnExp", value="logical"),
                 function(object, value){
                     if (length(value) == 1)
                         value <- rep(value, length(object))
                     if (length(object) != length(value))
                         stop("Length of replacement value is different than number of spectra.")
                     fData(object)$centroided <- value
                     if (validObject(object))
                         return(object)
                 })

############################################################
## rtime
##
## Get the retention time
setMethod("rtime", "OnDiskMSnExp",
          function(object){
              vals <- fData(object)$retentionTime
              names(vals) <- featureNames(object)
              return(vals)
          })

############################################################
## tic
##
## Get the total ion Current (from the featureData, thus it will not
## be recalculated).
setMethod("tic", "OnDiskMSnExp",
          function(object){
              vals <- fData(object)$totIonCurrent
              names(vals) <- featureNames(object)
              return(vals)
          })

############################################################
## polarity
##
## Get the polarity.
setMethod("polarity", "OnDiskMSnExp",
          function(object){
              vals <- fData(object)$polarity
              names(vals) <- featureNames(object)
              return(vals)
          })

#############################################################
## peaksCount
##
## Extract the peaksCount data; if the spectraProcessingQueue is empty
## we're just returning the peaksCount from the featureData, otherwise we
## check the function in there and eventually load the data and apply it.
setMethod("peaksCount",
          signature(object = "OnDiskMSnExp", scans = "missing"),
          function(object, scans, BPPARAM = bpparam()) {
              ## The feature data contains the original peaks
              ## count. This method fetches the peaks count from the
              ## (possibly processed) spectra.
              ##
              ## Add methods here that would not require raw data reading.
              skipFun <- c("removePeaks", "clean", "pickPeaks")


              if (length(object@spectraProcessingQueue) > 0) {
                  recalc <- any(unlist(lapply(object@spectraProcessingQueue,
                                              function(z) {
                                                  if (any(z@FUN == skipFun))
                                                      return(FALSE)
                                                  return(TRUE)
                                              })))
              } else {
                  ## No need to calculate the peaks count; we can use the
                  ## information from the feature data.
                  recalc <- FALSE
              }
              ## An important point here is that we DON'T want to get
              ## all of the data from all files in one go; that would
              ## require eventually lots of memory! It's better to do
              ## that per file; that way we could also do that in
              ## parallel.
              if (recalc) {
                  vals <- unlist(spectrapply(object,
                                  FUN = peaksCount,
                                  index = numeric(),
                                  BPPARAM = BPPARAM))
              } else {
                  vals <- fData(object)$originalPeaksCount
              }
              names(vals) <- featureNames(object)
              return(vals)
          })

############################################################
## ionCount
##
## Calculate the ion count, i.e. the sum of intensities per spectrum.
setMethod("ionCount", "OnDiskMSnExp",
          function(object, BPPARAM = bpparam()) {
              ## An important point here is that we DON'T want to get
              ## all of the data from all files in one go; that would
              ## require eventually lots of memory! It's better to do
              ## that per file; that way we could also do that in
              ## parallel.
              vals <- spectrapply(object,
                                  FUN = function(y) return(sum(y@intensity)),
                                  BPPARAM = BPPARAM)
              return(unlist(vals))
          })

############################################################
## spectra
##
## Extract the spectra of an experiment by reading the raw data from
## the original files, applying processing steps from the queue.
setMethod("spectra",
          "OnDiskMSnExp",
          function(object, BPPARAM = bpparam()){
              return(spectrapply(object, BPPARAM = BPPARAM))
          })

############################################################
## assayData
##
## Read the full data, put it into an environment and return that.
setMethod("assayData", "OnDiskMSnExp", function(object){
    message("Loading data from original files.")
    return(list2env(spectra(object)))
})

############################################################
## intensity
##
## Extract the intensity values of individual spectra. This means we
## have to read all of the data, create Spectrum objects, apply
## eventual processing steps and return the intensities.
setMethod("intensity", "OnDiskMSnExp", function(object, BPPARAM=bpparam()){
    message("Loading data to extract intensity values.")
    return(spectrapply(object, FUN=intensity, BPPARAM=BPPARAM))
})


############################################################
## mz
##
## Extract the mz values of individual spectra.
setMethod("mz", "OnDiskMSnExp", function(object, BPPARAM=bpparam()){
    message("Loading data to extract M/Z values.")
    return(spectrapply(object, FUN=mz, BPPARAM=BPPARAM))
})

############################################################
## [[
##
## Extract individual spectra (single ones).
setMethod("[[", "OnDiskMSnExp",
          function(x, i, j = "missing", drop = "missing") {
              if (length(i) != 1)
                  stop("subscript out of bounds")
              if (is.character(i))
                  i <- match(i, featureNames(x))
              if (is.na(i))
                  stop("subscript out of bounds")
              return(spectra(x[i])[[1]])
          })

############################################################
## [
##
## Subset by [
setMethod("[", signature(x = "OnDiskMSnExp",
                         i = "logicalOrNumeric",
                         j = "missing",
                         drop = "missing"),
          function(x, i, j, drop) {
              if (!(is.logical(i) | is.numeric(i)))
                  stop("Subsetting works only with numeric or logical!")
              if (is.logical(i)) {
                  if (length(i) != nrow(fData(x)))
                      stop("If 'i' is logical its length has to match the number of spectra!")
                  i <- which(i)
              }
              i <- sort(i)  ## Force sorted!
              ## Now subset the featureData. The function will
              ## complain if i is outside the range.
              x@featureData <- .subsetFeatureDataBy(featureData(x), index = i)
              return(x)
          })

############################################################
## spectrapply
##
## That's the main method to apply functions to the object's spectra, or
## to just return a list with the spectra, if FUN is empty.
## Arguments; all are related to the sub-setting of the featureData.
## index: subset by numeric index, logical or character. Character indices are "forwarded"
##  to argument "name".
## scanIdx: a numeric is expected, specifying the scan index. If a character vector
##  is provided, it is "forwarded" to argument "name".
## scanIdxCol: the column of the featureData containing the scan indices.
## name: a character vector, matching is performed using the row names of fd.
## rtlim: a numeric of length 2 specifying the retention time window from which spectra
##  should be extracted.
setMethod("spectrapply", "OnDiskMSnExp",
          function(object, FUN = NULL, index = NULL,
                   scanIdx = NULL, name = NULL,
                   rtlim = NULL,
                   BPPARAM = bpparam(), ...) {
    ## Check if we would do better with serial processing:
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)

    isOK <- .validateFeatureDataForOnDiskMSnExp(fd)
    if(!is.null(isOK))
        stop(isOK)
    fDataPerFile <- split(fd, f=fd$fileIdx)
    vals <- bplapply(fDataPerFile, FUN=.applyFun2SpectraOfFileMulti,
                     filenames=fileNames(object),
                     queue=object@spectraProcessingQueue,
                     APPLYFUN=FUN,
                     BPPARAM=BPPARAM)
    names(vals) <- NULL
    vals <- unlist(vals, recursive=FALSE)
              return(vals[rownames(fd)])
})

##============================================================
##  --  DATA MANIPULATION METHODS
##
##------------------------------------------------------------

############################################################
## removePeaks
##
## Add a "removePeaks" ProcessingStep to the queue and update
## the processingData information of the object.
setMethod("removePeaks", signature("OnDiskMSnExp"),
          function(object, t="min", verbose=TRUE){
              if(missing(t))
                  t <- "min"
              if(!is.numeric(t)){
                  if(t != "min")
                      stop("Argument 't' has to be either numeric or 'min'!")
              }
              ps <- ProcessingStep("removePeaks", list(t=t))
              ## Append the processing step to the queue.
              if(verbose)
                  message("Adding 'removePeaks' to the processing queue.")
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@removedPeaks <- c(object@processingData@removedPeaks,
                                                      as.character(t))
              object@processingData@processing <- c(object@processingData@processing,
                                                    paste("Curves <= ",
                                                          t
                                                         ," set to '0': ",
                                                          date(),sep=""))
              return(object)
          })

############################################################
## clean
##
## Add a "clean" ProcessingStep to the queue and update
## the processingData information of the object.
setMethod("clean", signature("OnDiskMSnExp"),
          function(object, all=FALSE, verbose=TRUE){
              if(!is.logical(all))
                  stop("Argument 'all' is supposed to be a logical!")
              ps <- ProcessingStep("clean", list(all=all))
              message("Adding 'clean' to the processing queue.")
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@cleaned <- TRUE
              object@processingData@processing <- c(object@processingData@processing,
                                                    paste0("Spectra cleaned: ", date()))
              return(object)
          })

############################################################
## trimMz
##
## Add the "trimMz" ProcessingStep to the queue and update the
## processingData.
setMethod("trimMz", signature("OnDiskMSnExp", "numeric"),
          function(object, mzlim, ...){
              ## Simple check on mzlim.
              if(length(mzlim) != 2)
                  stop("Argument 'mzlim' should be a numeric vector of length 2",
                       " specifying the lower and upper M/Z value range!")
              ps <- ProcessingStep("trimMz", list(mzlim=mzlim))
              message("Adding 'trimMz' to the processing queue.")
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              trmd <- object@processingData@trimmed
              ifelse(length(trmd)==0,
                     object@processingData@trimmed <- mzlim,
                     object@processingData@trimmed <- c(max(trmd[1],mzlim[1]),
                                                        min(trmd[2],mzlim[2])))
              object@processingData@processing <- c(object@processingData@processing,
                                                    paste("MZ trimmed [",
                                                          object@processingData@trimmed[1],
                                                          "..",
                                                          object@processingData@trimmed[2],
                                                          "]",sep=""))
              return(object)
          })

############################################################
## normalize
##
## Handle the 'normalize' method for MSnExp objects (calls normalise_MSnExp, and applies
## the normalization to each spectrum separately). Again we're adding a ProcessingStep
## for later, lazy processing.
setMethod("normalize", "OnDiskMSnExp",
          function(object, method=c("max", "sum"), ...){
              method <- match.arg(method)
              ps <- ProcessingStep("normalise", list(method=method))
              message("Adding 'normalize' to the processing queue.")
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@processing <- c(object@processingData@processing,
                                                    paste0("Spectra normalised (",method,"): ",
                                                           date()))
              object@processingData@normalised <- TRUE
              return(object)
          })

##============================================================
##  --  HELPER FUNCTIONS  --
##
##------------------------------------------------------------

############################################################
## .subsetFeatureDataBy
##
## Convenience function to subset a OnDiskMSnExp featureData based on
## provided subsetting criteria.
## index: subset by numeric index, logical or character. Character
##        indices are "forwarded" to argument "name".
## scanIdx: a numeric is expected, specifying the scan index. If a
##          character vector is provided, it is "forwarded" to
##          argument "name".
## scanIdxCol: the column of the featureData containing the scan
##             indices.
## name: a character vector, matching is performed using the row names
##       of fd.
## rtlim: a numeric of length 2 specifying the retention time window
##        from which spectra should be extracted.
.subsetFeatureDataBy <- function(fd, index=NULL, scanIdx=NULL, scanIdxCol="acquisitionNum",
                                 name=NULL, rtlim=NULL){
    ## First check index.
    if (length(index) > 0) {
        if (is.logical(index)) {
            if (length(index) != nrow(fd))
                stop("If 'index' is a logical vector its length has to match the number of",
                     " rows of the featureData!")
            index <- which(index)
        }
        if (is.numeric(index)) {
            gotIt <- index %in% 1:nrow(fd)
            if (!any(gotIt))
                stop("Provided indices are outside of the allowed range.")
            if (any(!gotIt))
                warning("Some of the provided indices are outside of the allowed range.")
            index <- index[gotIt]
            return(fd[index, , drop = FALSE])
        }
        if (is.character(index))
            name <- index
    }
    ## scanIdx
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
    ## name: subset by name, match to rownames.
    if (length(name) > 0) {
        gotIt <- name %in% rownames(fd)
        if (!any(gotIt))
            stop("None of the provided names found.")
        if (!all(gotIt))
            warning("Some of the provided names do not match featureData rownames.")
        name <- name[gotIt]
        return(fd[name, , drop=FALSE])
    }
    ## rtlim: subset by retention time range.
    if (length(rtlim > 0)){
        if (length(rtlim) > 2 | !is.numeric(rtlim))
            stop("Argument 'rtlim' has to be a numeric vector of length 2 specifying",
                 " the retention time window (range).")
        gotIt <- which(fd$retentionTime >= rtlim[1] & fd$retentionTime <= rtlim[2])
        if (length(gotIt) == 0)
            stop("No spectrum within the specified retention time window.")
        fd <- fd[gotIt, , drop=FALSE]
    }
    return(fd)
}

## Apply a function to the spectra of a file.
## The function will first read the raw data, create Spectrum objects from it, apply all
## ProcessingSteps apply the specified function and return its results. This is done on a
## per-file basis, so we can do that in parallel.
## Input args:
## fData: the data.frame containing all feature data info for the current file. All spectra
##   specified by column "spIdx" are read.
## filenames: the files from which the data should be read; this should be all file names
##   (i.e. those returned) by fileNames(object), as we're using the fileIdx column in fData
##   to get the "real" filename.
## queue: a list of ProcessingStep objects that should be applied to the Spectrum.
## APPLYFUN: the function that should be applied; if APPLYFUN is NULL it returns the spectra!
.applyFun2SpectraOfFile <- function(fData, filenames, queue=NULL,
                                    APPLYFUN=NULL){
    if(missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are requierd!")
    ## Subset the fData based on scnIDx; if provided.
    filename <- filenames[fData[1, "fileIdx"]]
    if(any(fData$msLevel > 1))
        stop("on-the-fly import currently only supported for MS1 level data.")
    ## Open the file.
    message("Reading data from file ", basename(filename), "...", appendLF=FALSE)
    fileh <- mzR::openMSfile(filename)
    on.exit(expr=mzR::close(fileh))
    on.exit(expr=message("OK"), add=TRUE)
    ## Now run the stuff per spectrum, i.e. read data, create object, apply the fun.
    ## Note that we're splitting the matrix not the data.frame, as that's faster.
    res <- lapply(split(fData[, c("fileIdx", "spIdx", "centroided", "acquisitionNum",
                                  "peaksCount", "totIonCurrent", "retentionTime",
                                  "polarity")], rownames(fData)),
                  FUN=function(z, theQ, fh, APPLYFUN){
                      ## Read the data.
                      spD <- mzR::peaks(fh, z[1, 2])
                      sp <- Spectrum1(peaksCount=z[1, 5], rt=z[1, 7],
                                      acquisitionNum=z[1, 4], scanIndex=z[1, 2],
                                      mz=spD[, 1], intensity=spD[, 2],
                                      fromFile=z[1, 1], centroided=z[1, 3],
                                      polarity=z[1, 8])
                      ## Now, apply the Queue.
                      if(length(theQ) > 0){
                          for(pStep in theQ){
                              sp <- execute(pStep, sp)
                          }
                      }
                      if(is.null(APPLYFUN))
                          return(sp)
                      ## Now what remains is to apply the APPLYFUN and return results.
                      return(APPLYFUN(sp))
                  }, theQ=queue, fh=fileh, APPLYFUN=APPLYFUN)
    return(res)
}
## The same function but using the standard "new" constructor from R.
.applyFun2SpectraOfFileSlow <- function(fData, filenames, queue=NULL,
                                        APPLYFUN=NULL){
    if(missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are requierd!")
    filename <- filenames[fData[1, "fileIdx"]]
    if(any(fData$msLevel > 1))
        stop("on-the-fly import currently only supported for MS1 level data.")
    ## Open the file.
    message("Reading data from file ", basename(filename), "...", appendLF=FALSE)
    fileh <- mzR::openMSfile(filename)
    on.exit(expr=mzR::close(fileh))
    on.exit(expr=message("OK"), add=TRUE)
    ## Now run the stuff per spectrum, i.e. read data, create object, apply the fun.
    ## Note that we're splitting the matrix not the data.frame, as that's faster.
    res <- lapply(split(fData[, c("fileIdx", "spIdx", "centroided", "acquisitionNum",
                                  "peaksCount", "totIonCurrent", "retentionTime",
                                  "polarity")], rownames(fData)),
                  FUN=function(z, theQ, fh, APPLYFUN){
                      ## Read the data.
                      spD <- mzR::peaks(fh, z[1, 2])
                      sp <- new("Spectrum1",
                                fromFile=z[1,1],
                                scanIndex=z[1, 2],
                                centroided=z[1,3],
                                acquisitionNum=z[1,4],
                                peaksCount=z[1,5],
                                rt=z[1,7],
                                polarity=z[1,8],
                                mz=spD[, 1],
                                intensity=spD[, 2])
                      ## Now, apply the Queue.
                      if(length(theQ) > 0){
                          for(i in 1:length(theQ)){
                              sp <- execute(theQ[[i]], sp)
                          }
                      }
                      if(is.null(APPLYFUN))
                          return(sp)
                      ## Now what remains is to apply the APPLYFUN and return results.
                      return(APPLYFUN(sp))
                  }, theQ=queue, fh=fileh, APPLYFUN=APPLYFUN)
    return(res)
}
## Using the C constructor that takes all values at once and creates a list of
## Spectrum1 objects, applies processing steps, applies the provided function and returns its
## results - or the list of Spectrum1 objects if APPLYFUN=NULL.
## Arguments:
## o fData: either a full data.frame (returned by fData(OnDiskMSnExp)) or a sub-set for
##   specific spectra. The data.frame should ONLY CONTAIN VALUES FOR SPECTRA OF ONE FILE!
## o filenames: fileNames(object)
## o queue: object@spectraProcessingQueue; if lenght > 0 all processing steps will be
##   applied to the created Spectrum1 objects.
## o APPLYFUN: the function to be applied to the Spectrum1 objects (such as ionCount etc).
##   If NULL the function returns the list of Spectrum1 objects.
.applyFun2SpectraOfFileMulti <- function(fData, filenames, queue=NULL,
                                         APPLYFUN=NULL){
    if(missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are required!")
    filename <- filenames[fData[1, "fileIdx"]]
    ## if(any(fData$msLevel > 1))
    ##     stop("on-the-fly import currently only supported for MS1 level data.")
    ## Open the file.
    message("Read data from file ", basename(filename), ".")
    fileh <- mzR::openMSfile(filename)
    on.exit(expr=mzR::close(fileh))
    msLevel1 <- which(fData$msLevel == 1)
    msLevelN <- which(fData$msLevel > 1)
    ## Process MS1 and MSn separately
    if (length(msLevel1) >= 1) {
        ms1fd <- fData[msLevel1, , drop = FALSE]
        ## Reading all of the data in "one go".
        allSpect <- mzR::peaks(fileh, ms1fd$spIdx)
        ## If we have more than one spectrum the peaks function returns a list.
        if (is(allSpect, "list")) {
            nValues <- lengths(allSpect) / 2
            allSpect <- do.call(rbind, allSpect)
        } else {
            ## otherwise it's a matrix, e.g. if only a single scan index was provided.
            nValues <- nrow(allSpect)
        }
        ## Call the C-constructor to create a list of Spectrum1 objects.
        res <- Spectra1(peaksCount = nValues,
                        scanIndex = ms1fd$spIdx,
                        rt = ms1fd$retentionTime,
                        acquisitionNum = ms1fd$acquisitionNum,
                        mz = allSpect[, 1],
                        intensity = allSpect[, 2],
                        centroided = ms1fd$centroided,
                        fromFile = ms1fd$fileIdx,
                        polarity = ms1fd$polarity,
                        nvalues = nValues)
        names(res) <- rownames(ms1fd)
    } else {
        res <- list()
    }
    if (length(msLevelN) >= 1) {
        msnfd <- split(fData[msLevelN, , drop = FALSE], f = 1:length(msLevelN))
        ## TODO: write/use C-constructor
        ## For now we're using the lapply, new() approach iteratively reading each
        ## spectrum from file and creating the Spectrum2.
        res2 <- lapply(msnfd, function(z) {
            spectD <- mzR::peaks(fileh, z$spIdx)
            return(new("Spectrum2",
                       merged = z$mergedScan,
                       precScanNum = z$precursorScanNum,
                       precursorMz = z$precursorMZ,
                       precursorIntensity = z$precursorIntensity,
                       precursorCharge = z$precursorCharge,
                       collisionEnergy = z$collisionEnergy,
                       msLevel = z$msLevel,
                       rt = z$retentionTime,
                       acquisitionNum = z$acquisitionNum,
                       scanIndex = z$spIdx,
                       mz = spectD[, 1],
                       intensity = spectD[, 2],
                       fromFile = z$fileIdx,
                       centroided = z$centroided,
                       polarity = z$polarity)
                   )
        })
        names(res2) <- rownames(fData)[msLevelN]
        res <- c(res, res2)
    }
    ## Ensure that ordering is the same than in fData:
    res <- res[match(rownames(fData), names(res))]
    ## If we have a non-empty queue, we might want to execute that too.
    if(!is.null(APPLYFUN) | length(queue) > 0){
        if(length(queue) > 0){
            message("Apply lazy processing steps:")
            for(j in 1:length(queue))
                message(" o '", queue[[j]]@FUN, "' with ", length(queue[[j]]@ARGS), " argument(s).")
        }
        res <- lapply(res, function(z, theQ, APPLF){
            if(length(theQ) > 0){
                for(pStep in theQ){
                    z <- execute(pStep, z)
                }
            }
            if(is.null(APPLF))
                return(z)
            return(APPLF(z))
        }, theQ=queue, APPLF=APPLYFUN)
    }
    message("DONE (", basename(filename), ")")
    return(res)
}

## Returns either NULL or a character string.
.validateFeatureDataForOnDiskMSnExp <- function(x){
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
