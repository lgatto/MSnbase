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

############################################################
## show
##
setMethod("show", "OnDiskMSnExp", function(object){
    procs <- processingQueue(object)
    ## Gather information.
    msLevels <- unique(msLevel(object))
    msnRt <- unname(rtime(object))
    nrt <- length(msnRt)
    rtr <- range(msnRt)
    cat("Object of class \"",class(object),"\"\n",sep="")
    cat(" Object size in memory: ")
    sz <- object.size(object)
    cat(round(sz/(1024^2),2),"Mb\n")
    cat("- - - Spectra data - - -\n")
    if (nrow(fData(object)) == 0) {
        cat(" none\n")
    } else {
        cat(" MS level(s):",msLevels,"\n")
        if (all(msLevel(object) > 1)) {
            ## cat(" Number of MS1 acquisitions:",nPrecScans,"\n")
            ## cat(" Number of MSn scans:",length(ls(assayData(object))),"\n")
            ## cat(" Number of precursor ions:",nPrecMz,"\n")
            ## if (nPrecMz > 0) {
            ##     cat("",uPrecMz,"unique MZs\n")
            ##     cat(" Precursor MZ's:",paste(signif(rangePrecMz,5),collapse=" - "),"\n")
            ## }
            ## cat(" MSn M/Z range:",round(msnMzRange,2),"\n")
        } else {
            cat(" Number of MS1 scans:",length(msnRt),"\n")
        }
        if (nrt > 0) {
            cat(" MSn retention times:",formatRt(rtr[1]),"-",formatRt(rtr[2]),"minutes\n")
        }
    }
    show(processingData(object))
    cat("- - - Meta data  - - -\n")
    Biobase:::.showAnnotatedDataFrame(phenoData(object),
                                      labels=list(object="phenoData"))
    cat("Loaded from:\n")
    f <- basename(processingData(object)@files)
    nf <- length(f)
    if (nf > 0) {
        if (nf < 3) {
            cat(paste0("  ", f, collapse = ", "), "\n")
        } else {
            cat("  [1]", paste(f[1], collapse = ", "))
            cat("...")
            cat("  [", nf, "] ", paste(f[nf], collapse = ", "),
                "\n", sep = "")
            cat("  Use 'fileNames(.)' to see all files.\n")
        }
    } else {
        cat(" none\n")
    }
    Biobase:::.showAnnotatedDataFrame(protocolData(object),
                                      labels=list(object="protocolData"))
    Biobase:::.showAnnotatedDataFrame(featureData(object),
                                      labels=list(
                                          object="featureData",
                                          sampleNames="featureNames",
                                          varLabels="fvarLabels",
                                          varMetadata="fvarMetadata"))
    cat("experimentData: use 'experimentData(object)'\n")
    pmids <- pubMedIds(object)
    if (length(pmids) > 0 && all(pmids != ""))
        cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
    invisible(NULL)
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
## header
##
## Extract the "header" data from the featureData, rename some
## of the columns and return the data.frame.
setMethod("header",
          signature("OnDiskMSnExp","missing"),
          function(object){
              hd <- fData(object)
              ## If we need the "ionCount" column we might have to read it on-the-fly,
              ## which would badly suck, as it would slow down stuff.
              ## Rename columns
              colnames(hd)[colnames(hd) == "fileIdx"] <- "file"
              colnames(hd)[colnames(hd) == "retentionTime"] <- "retention.time"
              colnames(hd)[colnames(hd) == "peaksCount"] <- "peaks.count"
              colnames(hd)[colnames(hd) == "totIonCurrent"] <- "tic"
              colnames(hd)[colnames(hd) == "msLevel"] <- "ms.level"
              colnames(hd)[colnames(hd) == "acquisitionNum"] <- "acquisition.number"
              return(hd)
          })
setMethod("header",
          signature=c("OnDiskMSnExp","numeric"),
          function(object, scans){
              hd <- .subsetFeatureDataBy(fData(object), scanIdx=scans)
              colnames(hd)[colnames(hd) == "fileIdx"] <- "file"
              colnames(hd)[colnames(hd) == "retentionTime"] <- "retention.time"
              colnames(hd)[colnames(hd) == "peaksCount"] <- "peaks.count"
              colnames(hd)[colnames(hd) == "totIonCurrent"] <- "tic"
              colnames(hd)[colnames(hd) == "msLevel"] <- "ms.level"
              colnames(hd)[colnames(hd) == "acquisitionNum"] <- "acquisition.number"
              return(hd)
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
setMethod("peaksCount", signature(object="OnDiskMSnExp", scans="missing"),
          function(object, scans, BPPARAM=bpparam()){
              scans <- numeric()
              return(peaksCount(object, scans=scans, BPPARAM=BPPARAM))
          })
setMethod("peaksCount", signature(object="OnDiskMSnExp", scans="numeric"),
          function(object, scans, BPPARAM=bpparam()){
              fd <- .subsetFeatureDataBy(fData(object), index=scans)
              ## Peaks count from the original files is available in the featureData,
              ## thus, we only have to read the raw data again and calculate the peaksCount
              ## if we did do some processing of the data. And here also just processing
              ## that has an influence on the number of data duplets. So, "removePeaks" does
              ## not change the peaksCount, only "clean" would. We thus check the methods
              ## in 'spectraProcessingQueue' and decide whether or not we have to load
              ## the data.
              skipFun <- c("removePeaks")  ## Add methods here that would not require raw data reading.
              recalc <- FALSE
              if(length(object@spectraProcessingQueue) > 0){
                  recalc <- any(unlist(lapply(object@spectraProcessingQueue,
                                              function(z){
                                                  if(any(z@FUN == skipFun))
                                                      return(FALSE)
                                                  return(TRUE)
                                              })))
              }
              if(recalc){
                  ## Heck; reload the data.
                  message("Loading the raw data to calculate peaksCount.")
                  ## An important point here is that we DON'T want to get all of the data from
                  ## all files in one go; that would require eventually lots of memory! It's better
                  ## to do that per file; that way we could also do that in parallel.
                  fDataPerFile <- split(fd, f=fd$fileIdx)
                  vals <- bplapply(fDataPerFile, FUN=.applyFun2SpectraOfFileMulti,
                                   filenames=fileNames(object),
                                   queue=object@spectraProcessingQueue,
                                   APPLYFUN=peaksCount,
                                   BPPARAM=BPPARAM)
                  names(vals) <- NULL
                  vals <- unlist(vals, recursive=TRUE)
                  return(vals[rownames(fd)])
              }else{
                  vals <- fd$peaksCount
                  names(vals) <- rownames(fd)
              }
              return(vals)
          })

############################################################
## ionCount
##
## Calculate the ion count, i.e. the sum of intensities per spectrum.
setMethod("ionCount", "OnDiskMSnExp",
          function(object, BPPARAM=bpparam()){
              fd <- fData(object)
              ## Heck; reload the data.
              message("Loading the raw data to calculate peaksCount.")
              ## An important point here is that we DON'T want to get all of the data from
              ## all files in one go; that would require eventually lots of memory! It's better
              ## to do that per file; that way we could also do that in parallel.
              fDataPerFile <- split(fd, f=fd$fileIdx)
              vals <- bplapply(fDataPerFile, FUN=.applyFun2SpectraOfFileMulti,
                               filenames=fileNames(object),
                               queue=object@spectraProcessingQueue,
                               APPLYFUN=function(y){return(sum(y@intensity))},
                               BPPARAM=BPPARAM)
              names(vals) <- NULL
              vals <- unlist(vals, recursive=TRUE)
              return(vals[rownames(fd)])
          })

############################################################
## spectra
##
## Extract the spectra of an experiment by reading the raw data from
## the original files, applying processing steps from the queue.
## The optional argument "scans" allows to restrict extraction of
## specific scans.
## Returned spectra are always sorted by scan index and file (first scan
## first file, first scan second file, etc.).
setMethod("spectra", "OnDiskMSnExp", function(object, scans=NULL,
                                              BPPARAM=bpparam()){
    fd <- .subsetFeatureDataBy(fData(object), index=scans)
    fdPerFile <- split(fd, f=fd$fileIdx)
    ## Overwriting the BPPARAM setting if we've only got some scans!
    if(length(scans) > 0 & length(scans) < 800)
        BPPARAM <- SerialParam()
    vals <- bplapply(fdPerFile, FUN=.applyFun2SpectraOfFileMulti,
                     filenames=fileNames(object),
                     queue=object@spectraProcessingQueue,
                     BPPARAM=BPPARAM)
    names(vals) <- NULL
    vals <- unlist(vals)
    return(vals[rownames(fd)])
})

############################################################
## assayData
##
## Read the full data, put it into an environment and return that.
setMethod("assayData", "OnDiskMSnExp", function(object){
    fd <- fData(object)
    fdPerFile <- split(fd, f=fd$fileIdx)
    vals <- bplapply(fdPerFile, FUN=.applyFun2SpectraOfFileMulti,
                     filenames=fileNames(object),
                     queue=object@spectraProcessingQueue,
                     BPPARAM=bpparam())
    names(vals) <- NULL
    vals <- unlist(vals)
    return(list2env(vals[rownames(fd)]))
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
              if(verbose)
                  message("Adding 'clean' to the processing queue.")
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@cleaned <- TRUE
              object@processingData@processing <- c(object@processingData@processing,
                                                    paste0("Spectra cleaned: ", date()))
              return(object)
          })


##============================================================
##  --  HELPER FUNCTIONS  --
##
##------------------------------------------------------------

############################################################
## .subsetFeatureDataBy
##
## Convenience function to subset a OnDiskMSnExp featureData based on provided
## subsetting criteria.
## index: subset by numeric index, logical or character. Character indices are "forwarded"
##  to argument "name".
## scanIdx: a numeric is expected, specifying the scan index. If a character vector
##  is provided, it is "forwarded" to argument "name".
## scanIdxCol: the column of the featureData containing the scan indices.
## name: a character vector, matching is performed using the row names of fd.
.subsetFeatureDataBy <- function(fd, index=NULL, scanIdx=NULL, scanIdxCol="spIdx",
                                 name=NULL){
    ## First check index.
    if(length(index) > 0){
        if(is.logical(index)){
            if(length(index) != nrow(fd))
                stop("If 'index' is a logical vector its length has to match the number of",
                     " rows of the featureData!")
            index <- which(index)
        }
        if(is.numeric(index)){
            gotIt <- index %in% 1:nrow(fd)
            if(!any(gotIt))
                stop("Provided indices are outside of the allowed range.")
            if(any(!gotIt))
                warning("Some of the provided indices are outside of the allowed range.")
            index <- index[gotIt]
            return(fd[index, , drop=FALSE])
        }
        if(is.character(index))
            name <- index
    }
    ## scanIdx
    if(length(scanIdx) > 0){
        if(is.numeric(scanIdx)){
            gotIt <- scanIdx %in% fd[, scanIdxCol]
            if(!any(gotIt))
                stop("None of the provided scan indices are available!")
            if(!all(gotIt))
                warning("Some of the provided scan indices are not available.")
            return(fd[which(fd[, scanIdxCol] %in% scanIdx), , drop=FALSE])
        }
        if(is.character(scanIdx))
            name <- scanIdx
    }
    ## name: subset by name, match to rownames.
    if(length(name) > 0){
        gotIt <- name %in% rownames(fd)
        if(!any(gotIt))
            stop("None of the provided names found.")
        if(!all(gotIt))
            warning("Some of the provided names do not match featureData rownames.")
        name <- name[gotIt]
        return(fd[name, , drop=FALSE])
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
                  FUN=function(z, theQ, skipValidate, fh, APPLYFUN){
                      ## Read the data.
                      spD <- mzR::peaks(fh, z[1, 2])
                      sp <- Spectrum1(peaksCount=z[1, 5], rt=z[1, 7],
                                      acquisitionNum=z[1, 4], scanIndex=z[1, 2],
                                      tic=z[1, 6], mz=spD[, 1], intensity=spD[, 2],
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
                  }, theQ=queue, skipValidate=fast, fh=fileh, APPLYFUN=APPLYFUN)
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
                  FUN=function(z, theQ, skipValidate, fh, APPLYFUN){
                      ## Read the data.
                      spD <- mzR::peaks(fh, z[1, 2])
                      sp <- new("Spectrum1",
                                fromFile=z[1,1],
                                scanIndex=z[1, 2],
                                centroided=z[1,3],
                                acquisitionNum=z[1,4],
                                peaksCount=z[1,5],
                                tic=z[1,6],
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
                  }, theQ=queue, skipValidate=fast, fh=fileh, APPLYFUN=APPLYFUN)
    return(res)
}
## Using the C constructor that takes all values and creates a list of Spectrum1 objects.
.applyFun2SpectraOfFileMulti <- function(fData, filenames, queue=NULL,
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
    ## Reading all of the data in "one go".
    allSpect <- mzR::peaks(fileh, fData$spIdx)
    ## If we have more than one spectrum.
    if(is(allSpect, "list")){
        nValues <- lengths(allSpect) / 2
        allSpect <- do.call(rbind, allSpect)
    }else{
        ## otherwise it's a matrix, e.g. if only a single scan index was provided.
        nValues <- nrow(allSpect)
    }
    ## Call the C-constructor to create a list of Spectrum1 objects.
    res <- Spectra1(peaksCount=fData$peaksCount,
                    scanIndex=fData$spIdx,
                    rt=fData$retentionTime,
                    acquisitionNum=fData$acquisitionNum,
                    tic=fData$totIonCurrent,
                    mz=allSpect[, 1],
                    intensity=allSpect[, 2],
                    centroided=fData$centroided,
                    fromFile=fData$fileIdx,
                    polarity=fData$polarity,
                    nvalues=nValues)
    names(res) <- rownames(fData)

    ## If we have a non-empty queue, we might want to execute that too.
    if(!is.null(APPLYFUN) | length(queue) > 0){
        cat("APPLYFUN not empty\n")
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
    message("OK")
    return(res)
}

