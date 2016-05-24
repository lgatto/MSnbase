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
## Extract the peaksCount data; if the spectraProcessingQueue is empty
## we're just returning the peaksCount from the featureData, otherwise we
## check the function in there and eventually load the data and apply it.
setMethod("peaksCount", signature(object="OnDiskMSnExp", scans="missing"),
          function(object, scans, method=1, BPPARAM=bpparam()){
              ## Now, that's a little different now; we might have to calculate that
              ## if we've done some data processing. Thus we have to check the
              ## spectraProcessingQueue.
              ## By default we will always re-calculate the peaks count, but we define
              ## some functions from which we know that they don't change that count.
              ## removePeaks for example does NOT remove peaks, onyl "clean" does. So, if
              ## there are some new methods that do NOT remove peaks they should be added
              ## to this list.
              skipFun <- c("removePeaks")
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
                  ## peaks count is the length of mz values per spectrum.
                  ## The "clean", but somewhat slow version is to get all of the spectra,
                  ## run the processing function on these spectra and call peaksCount on them.
                  fDataPerFile <- split(fData(object), f=fData(object)$fileIdx)
                  if(method == 1){
                      message("Apply on Spectrum")
                      vals <- lapply(fDataPerFile, FUN=.applyFun2SpectraOfFile,
                                     filenames=fileNames(object),
                                     queue=object@spectraProcessingQueue,
                                     APPLYFUN=peaksCount)
                  }
                  if(method == 2){
                      message("Using the multi-constructor in C")
                      vals <- lapply(fDataPerFile, FUN=.applyFun2SpectraOfFileMulti,
                                     filenames=fileNames(object),
                                     queue=object@spectraProcessingQueue,
                                     APPLYFUN=peaksCount)
                  }
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
## Extract the spectra of an experiment by reading the raw data from
## the original files, applying processing steps from the queue.
## The optional argument "scans" allows to restrict extraction of
## specific scans.
## Returned spectra are always sorted by scan index and file (first scan
## first file, first scan second file, etc.).
setMethod("spectra", "OnDiskMSnExp", function(object, scans=NULL,
                                              method=1,
                                              BPPARAM=bpparam()){
    fd <- .subsetFeatureDataBy(fData(object), scanIdx=scans)
    ## fd <- fData(object)
    ## ## Check scans argument.
    ## if(length(scans) > 0){
    ##     if(!is.numeric(scans))
    ##         stop("Argument 'scans', if provided, has to be numeric!")
    ##     ## Subset the featureData based on scans.
    ##     fd <- fd[fd$spIdx %in% scans, , drop=FALSE]
    ##     if(nrow(fd) == 0)
    ##         stop("The requested scan indices are not available!")
    ## }
    fdPerFile <- split(fd, f=fd$fileIdx)
    if(method == 1){
        vals <- lapply(fdPerFile, FUN=.applyFun2SpectraOfFile,
                       filenames=fileNames(object),
                       queue=object@spectraProcessingQueue)
    }
    if(method == 2){
        vals <- lapply(fdPerFile, FUN=.applyFun2SpectraOfFileMulti,
                       filenames=fileNames(object),
                       queue=object@spectraProcessingQueue)
    }
    names(vals) <- NULL
    vals <- unlist(vals)
    return(vals[rownames(fd)])
})

############################################################
##  --  HELPER FUNCTIONS  --
##
############################################################

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

############################################################
## .getDataDupletsAsMatrix
##
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
    res <- mapply(fileNames(object), wantScnIdx, FUN=.getDataDupletsAsMatrixForFile,
                  SIMPLIFY=FALSE)
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

    nValues <- lengths(allSpect) / 2
    allSpect <- do.call(rbind, allSpect)
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
        res <- lapply(res, function(z, theQ, APPLYFUN){
            if(length(theQ) > 0){
                for(pStep in theQ){
                    z <- execute(pStep, z)
                }
                if(is.null(APPLYFUN))
                    return(z)
                return(APPLYFUN(z))
            }
        }, theQ=queue, APPLYFUN=APPLYFUN)
    }
    message("OK")
    return(res)
}

