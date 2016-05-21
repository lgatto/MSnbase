############################################################
## Methods for OnDiskMSnExp objects.

setMethod("initialize",
          signature(.Object="OnDiskMSnExp"),
          function(.Object, ...){
              ## if(!any(names(list(...)) == "onDisk")){
              ##     .Object@onDisk=TRUE
              ## }
              callNextMethod()
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
              hd <- header(object)
              return(hd[scans, ])
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
## peaksCount
##
## Extract the peaksCount data.
setMethod("peaksCount", signature(object="OnDiskMSnExp", scans="missing"),
          function(object, scans){
              ## Now, that's a little different now; we might have to calculate that
              ## if we've done some data processing. Thus we have to check the
              ## spectraProcessingQueue.
              ## Define a list of function names that, if present, require a re-calculation
              ## of the peaksCount.
              ## Note: only if we have clean we will have to re-count peaks! removePeaks
              ## does not remove peaks, it does only flood them.
              considerFun <- c("removePeaks", "clean")
              recalc <- FALSE
              if(length(object@spectraProcessingQueue) > 0){
                  recalc <- any(unlist(lapply(object@spectraProcessingQueue,
                                               function(z){
                                                   return(any(z@FUN == considerFun))
                                               })))
              }
              if(recalc){
                  ## Heck; reload the data.
                  message("Loading the raw data to calculate peaksCount.")
                  ## peaks count is the length of mz values per spectrum.
                  ## The "clean", but somewhat slow version is to get all of the spectra,
                  ## run the processing function on these spectra and call peaksCount on them.
                  fDataPerFile <- split(fData(object), f=fData(object)$fileIdx)
                  vals <- lapply(fDataPerFile, FUN=.applyFun2SpectraOfFile, filenames=fileNames(object),
                                 queue=object@spectraProcessingQueue, APPLYFUN=peaksCount)
                  names(vals) <- NULL
                  vals <- unlist(vals, recursive=TRUE)
                  return(vals[featureNames(object)])
                  ## An important point here is that we DON'T want to get all of the data from
                  ## all files in one go; that would require eventually lots of memory! It's better
                  ## to do that per file; that way we could also do that in parallel.
              }else{
                  vals <- fData(object)$peaksCount
                  names(vals) <- featureNames(object)
              }
              return(vals)
          })

############################################################
## spectra
##
## Extract the spectra.
setMethod("spectra", "OnDiskMSnExp", function(object, ...){
    fDataPerFile <- split(fData(object), f=fData(object)$fileIdx)
    vals <- lapply(fDataPerFile, FUN=.applyFun2SpectraOfFile, filenames=fileNames(object),
                   queue=object@spectraProcessingQueue)
    names(vals) <- NULL
    vals <- unlist(vals)
    return(vals[featureNames(object)])
})

############################################################
##  --  HELPER FUNCTIONS  --
##
############################################################

## Get the mz and intensity values as a matrix.
## object: the OnDiskMSnExp object.
## scnIdx: optionally the index of the scans that should be extracted. The function checks if these
##         indices are available and throws a warning (or error) if some or none are present.
##
## return: a list of data matrices. Names of the list elements will be the spectrum names.
.getDataDupletsAsMatrix <- function(object, scnIdx){
    ## fileNames(object) to get the file names; open the connection.
    ## scanIndex(object) to get all of the scan indices.
    ## fromFile(object) to get to know from which file the spectra should be extracted.
    ## Do that with PARALLEL processing.
    ## Split the scanIndex by from File.
    wantScnIdx <- scanIndex(object)
    wantScnIdx <- split(wantScnIdx, fromFile(object))
    if(!missing(scnIdx)){
        ## Subset the wantScnIdx and check if the requested scnIdx are available.
        wantScnIdx <- lapply(wantScnIdx, FUN=function(z){
            newz <- z[z %in% scnIdx]
            if(length(newz) == 0)
                stop("The requested scan indices are not available!")
            if(any(!(scnIdx %in% newz)))
                warning("Some of the requested scan indices are not available.")
            return(newz)
        })
    }
    ## Do the call:
    res <- mapply(fileNames(object), wantScnIdx, FUN=.getDataDupletsAsMatrixForFile, SIMPLIFY=FALSE)
    names(res) <- NULL
    res <- unlist(res, recursive=FALSE)
    return(res[featureNames(object)])
}
## The function we'll apply in the mapply (bplapply).
.getDataDupletsAsMatrixForFile <- function(filename, wantScnIdx){
    ## Open the handle.
    message("Reading data from file ", basename(filename), "...", appendLF=FALSE)
    fileh <- mzR::openMSfile(filename)
    ## Close handle
    on.exit(expr=mzR::close(fileh))
    on.exit(expr=message("OK"), add=TRUE)

    ## Extract the data.
    return(mzR::peaks(fileh, wantScnIdx))
}

## Apply a function to the spectra of a file.
## The function will first read the raw data, create Spectrum objects from it, apply all
## ProcessingSteps apply the specified function and return its results. This is done on a
## per-file basis, so we can do that in parallel.
## Input args:
## filenames: the files from which the data should be read; this should be all file names (i.e. those returned)
##            by fileNames(object), as we're using the fileIdx column in fData to get the "real" filename.
## fData: the data.frame containing all feature data info for the current file.
## scnIdx: optional index of the scans that should be extracted. The function checks if these
##         are present (in fData) and throws a warning (or error) if not (error if none is present).
## queue: a list of ProcessingStep objects that should be applied to the Spectrum.
## APPLYFUN: the function that should be applied; if APPLYFUN is NULL it returns the spectra!
## fast: if that's TRUE we're not going to do a validity check!
## Some notes on performance: splitting the featureData as matrix would be faster, but then we have to
##   call as.integer etc on the individual fields. Also, reading the full data and running the lapply on
##   that doesn't improve speed.
.applyFun2SpectraOfFile <- function(fData, scnIdx, filenames, queue=NULL, APPLYFUN=NULL, fast=TRUE){
    ## Subset the fData based on scnIDx; if provided.
    if(!missing(scnIdx)){
        gotIt <- fData$spIdx %in% scnIdx
        if(!any(gotIt))
            stop("None of the requested scan indices is present!")
        fData <- fData[gotIt, , drop=FALSE]
        if(any(!(scnIdx %in% fData$spIdx)))
            warning("Some of the requested scan indices are not available.")
    }
    filename <- filenames[fData[1, "fileIdx"]]
    if(any(fData$msLevel > 1))
        stop("on-the-fly import currently only supported for MS1 level data.")
    ## Open the file.
    message("Reading data from file ", basename(filename), "...", appendLF=FALSE)
    fileh <- mzR::openMSfile(filename)
    on.exit(expr=mzR::close(fileh))
    on.exit(expr=message("OK"), add=TRUE)
    ## Now run the stuff per spectrum, i.e. read data, create object, apply the fun. Note that we're
    ## splitting the matrix not the data.frame, as that's faster.
    res <- lapply(split(fData[, c("fileIdx", "spIdx", "centroided", "acquisitionNum", "peaksCount",
                                  "totIonCurrent", "retentionTime", "polarity")], rownames(fData)),
                  FUN=function(z, theQ, skipValidate, fh, APPLYFUN){
                      ## Read the data.
                      spD <- mzR::peaks(fh, z[1, 2])
                      if(skipValidate){
                          sp <- new("Spectrum1")
                          sp@fromFile <- z[1,1]
                          sp@scanIndex <- z[1, 2]
                          sp@centroided <- z[1,3]
                          sp@acquisitionNum <- z[1,4]
                          sp@peaksCount <- z[1,5]
                          sp@tic <- z[1,6]
                          sp@rt <- z[1,7]
                          sp@polarity <- z[1,8]
                          sp@mz <- spD[, 1]
                          sp@intensity <- spD[, 2]
                      }else{
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
                      }
                      ## We could speed up things if we're calling that *before* creating
                      ## the spectrum!
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


.applyFun2SpectraOfFile2 <- function(fData, scnIdx, filenames, queue=NULL, APPLYFUN=NULL, fast=TRUE){
    ## Subset the fData based on scnIDx; if provided.
    if(!missing(scnIdx)){
        gotIt <- fData$spIdx %in% scnIdx
        if(!any(gotIt))
            stop("None of the requested scan indices is present!")
        fData <- fData[gotIt, , drop=FALSE]
        if(any(!(scnIdx %in% fData$spIdx)))
            warning("Some of the requested scan indices are not available.")
    }
    filename <- filenames[fData[1, "fileIdx"]]
    if(any(fData$msLevel > 1))
        stop("on-the-fly import currently only supported for MS1 level data.")
    ## Open the file.
    message("Reading data from file ", basename(filename), "...", appendLF=FALSE)
    fileh <- mzR::openMSfile(filename)
    on.exit(expr=mzR::close(fileh))
    on.exit(expr=message("OK"), add=TRUE)
    ## Now run the stuff per spectrum, i.e. read data, create object, apply the fun. Note that we're
    ## splitting the matrix not the data.frame, as that's faster.
    res <- lapply(split(fData[, c("fileIdx", "spIdx", "centroided", "acquisitionNum", "peaksCount",
                                  "totIonCurrent", "retentionTime", "polarity")], rownames(fData)),
                  FUN=function(z, theQ, skipValidate, fh, APPLYFUN){
                      ## Read the data.
                      spD <- mzR::peaks(fh, z[1, 2])
                      sp <- Spectrum1(msLevel=1, peaksCount=z[1, 5], rt=z[1, 7], acquisitionNum=z[1, 4],
                                      scanIndex=z[1, 2], tic=z[1, 6], mz=spD[, 1], intensity=spD[, 2],
                                      fromFile=z[1, 1], centroided=z[1, 3], polarity=z[1, 8])
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


## Test that.
## fname <- fileNames(odmse)[1]
## fd <- fData(odmse)
## fd <- fd[fd$fileIdx == 1, ]

## system.time(
##     Test <- MSnbase:::.applyFun2SpectraOfFile(fname, fd)
## )  ## 13.5 secs.
## system.time(
##     Test2 <- MSnbase:::.applyFun2SpectraOfFile(fname, fd, fast=FALSE)
## )  ## 20 secs.
## system.time(
##     Test3 <- MSnbase:::.applyFun2SpectraOfFile2(fname, fd)
## )  ## 14.4
## system.time(
##     Test4 <- MSnbase:::.applyFun2SpectraOfFile2(fname, fd, fast=FALSE)
## )  ## 18.3
## system.time(
##     Test5 <- MSnbase:::.applyFun2SpectraOfFile3(fname, fd)
## )  ## 14.1
## system.time(
##     Test6 <- MSnbase:::.applyFun2SpectraOfFile3(fname, fd, fast=FALSE)
## )  ## 18.0

## library(RUnit)
## checkEquals(Test, Test2)
## checkEquals(Test2, Test3)
## checkEquals(Test3, Test4)
## checkEquals(Test4, Test5)
## checkEquals(Test5, Test6)
