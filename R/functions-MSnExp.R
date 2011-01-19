"[.MSnExp" <- function(x,i) {
  if (max(i)>length(x) | min(i)<1)
    stop("subscript out of bonds")
  whichElements <- ls(assayData(x))[i]
  x@assayData <- list2env(mget(whichElements,assayData(x)))
  featureData(x) <- featureData(x)[i,]
  return(x)
}

header.MSnExp <- function(object) {
  if (any(msLevel(object)<2))
    stop("header() only works for MS levels > 1.")
  tbl <- table(fromFile(object))
  idx <- as.numeric(unlist(apply(tbl,1,function(x) 1:x)))
  return(data.frame(cbind(index=idx,
                          file=fromFile(object),
                          retention.time=rtime(object),
                          precursor.mz=precursorMz(object),
                          peaks.count=peaksCount(object),
                          tic=tic(object),
                          ms.level=msLevel(object),
                          charge=precursorCharge(object),
                          collision.energy=collisionEnergy(object))))
       }

## mergeSpectra <- function(object, ## MSnExp object
##                          fun=sum,
##                          verbose=TRUE) {
##   spectra <- spectra(object)
##   prec <- sapply(spectra,function(x) x@precursorMz)
##   uprec <- unique(prec)
##   luprec <- length(uprec)
##   if (verbose) 
##     cat(luprec,"unique precursors for",length(spectra),"spectra\n")
##   newSpectra <- vector("list",length=luprec)
##   if (verbose)
##     pb <- txtProgressBar(min = 0, max = luprec, style = 3)
##   for (i in 1:luprec) {
##     if (verbose)
##       setTxtProgressBar(pb, i)   
##     pi <- uprec[i]
##     sel <- prec %in% pi
##     l <- spectra[sel]
##     ## unique M/Z values
##     mzs <- unique(unlist(lapply(l,function(x) x@mz)))
##     ## vectors of intensitues (ints) and number of times
##     ## a given itensity is found (icounts)
##     ints <- icounts <- numeric(length(mzs))
##     names(icounts) <- mzs
##     allints <- unlist(lapply(l,function(x) x@intensity))
##     names(allints) <- unlist(lapply(l,function(x) x@mz))
##     for (j in 1:length(mzs)) {
##       k <- as.character(mzs[j])
##       ints[j] <- fun(allints[names(allints) %in% k])
##       icounts[j] <- sum(names(allints) %in% k)
##     }
##     o <- order(mzs)
##     newSpectra[[i]] <- new("Spectrum",
##                            merged=which(sel),
##                            rt=unique(sapply(l,function(x) x@rt)),
##                            msLevel=unique(sapply(l,function(x) x@msLevel)),
##                            precursorMz=pi,
##                            peaksCount=length(ints),
##                            intensity=ints[o],
##                            mz=mzs[o])
##     validObject(newSpectra[[i]])
##   }
##   if (verbose)
##     close(pb)
##   object@process@processing <- c(object@process@processing,
##                                  paste("Precursor ions with identical M/Z merged:",
##                                        date()))
##   object@process@merged <- TRUE
##   ## TODO update and set featureData
##   return(new("MSnExp",
##              assayData=list2env(newSpectra),
##              metadata=object@metadata,
##              process=object@process))
## }

extractPrecSpectra <- function(object,prec) {
  sel <- sapply(spectra(object),function(x) x@precursorMz %in% prec)
  n <- length(prec)
  object@processingData@processing <- c(object@processingData@processing,
                                        paste(n,"precursors extracted:",date()))
  ## TODO: need to update phenoData - check files
  return(new("MSnExp",
             spectra=list2env(object@spectra[sel]),
             processingData=object@processingData,
             assayData=object@assayData,
             phenoData=object@phenoData,
             featureData=object@featureData[sel,],
             experimentData=object@experimentData,
             protocolData=object@protocolData))
}


extractSpectra <- function(object,selected) {
  if (!is.logical(selected))
    stop("Please use logicals to extract a set of spectra.")
  object@spectra <- spectra(object)[selected]
  object@fromFile <- object@fromFile[selected]
  return(object)
}

removePeaks.MSnExp <- function(object,t="min",verbose=TRUE) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <-  llply(spectra(object),function(x) removePeaks(x,t),.progress=progress)
  object@assayData <- list2env(spectraList)
  object@processingData@removedPeaks <- c(object@processingData@removedPeaks,
                                          as.character(t))
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Curves <= ",
                                              t
                                              ," set to '0': ",
                                              date(),sep=""))
  return(object)
}



clean.MSnExp <- function(object,verbose=TRUE) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectra <- llply(spectra(object),function(x) clean(x),.progress=progress)
  object@assayData <- list2env(spectra)
  object@processingData@cleaned <- TRUE
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Spectra cleaned: ",date(),sep=""))
  return(object)
}

quantify.MSnExp <- function(object,reporters,method,verbose) {
  ## Display progress bar with eapply
  ## if (verbose) {
  ##   ._cnt <- 1
  ##   pb <- txtProgressBar(min = 0, max = length(object), style = 3)
  ##   ## Quantification -- creating exprs for assayData slot
  ##   peakData <- eapply(assayData(object),function(x) {
  ##     setTxtProgressBar(pb, ._cnt)
  ##     ._cnt <<- ._cnt+1
  ##     quantify(x,reporters,method)
  ##   })
  ##   close(pb)
  ##   rm(pb)
  ##   rm(._cnt)
  ## } else {
  ##   peakData <- eapply(assayData(object),quantify,reporters,method)
  ## }

  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <- spectra(object)
  ## Quantification -- creating exprs for assayData slot
  if (length(spectraList)==1) {
    peakData <- quantify(spectraList[[1]],reporters,method)
    .exprs <- t(peakData$peakQuant)
    .qual <- t(peakData$curveData)
  } else {
    peakData <- llply(spectraList,quantify,reporters,method,.progress=progress)
    .exprs <- do.call(rbind,sapply(peakData,"[","peakQuant"))
    .qual <- do.call(rbind,sapply(peakData,"[","curveStats"))
  }
  rownames(.exprs) <- sub(".peakQuant","",rownames(.exprs))
  rownames(.qual) <- sub(".curveStats","",rownames(.qual))
  
  ## Updating MSnprocess slot
  object@processingData@processing <- c(object@processingData@processing,
                                        paste(reporters@name,
                                              " quantification by ",method,
                                              ": ",date(),sep=""))
  object@processingData@centroided <- TRUE
                 
  ## Creating new MSnSet
  msnset <- new("MSnSet",
                qual=.qual,
                exprs=.exprs, 
                processingData=object@processingData,
                experimentData=experimentData(object),
                ## protocolData=protocolData(object) ## updated/added below
                ## phenoData=pdata,                  ## updated/added below
                ## featureData=featureData(object),  ## updated/added below
                annotation="No annotation")
  
  ## Updating featureData slot or creating one
  fd <- header(object)
  if (nrow(fData(object))>0) { 
    if (nrow(fData(object))==length(object)) {
      fd <- combine(fData(object),fd)
    } else {
      warning("Unexpected number of features in featureData slot. Dropping it.")
    }
  }
  ## featureData rows must be reordered to match assayData rows
  .featureData <- new("AnnotatedDataFrame",data=fd[rownames(.exprs),])
  featureData(msnset) <- .featureData

  ## Updating phenoData slot or creating one
  .phenoData <- new("AnnotatedDataFrame",
                    data=data.frame(mz=reporters@mz,
                      reporters=reporters@name,
                      row.names=reporters@reporterNames))
  if (nrow(pData(object))>0) { 
    if (nrow(pData(object))==length(reporters)) {
      .phenoData <- combine(phenoData(object),.phenoData)
    } else {
      ## Here, something more clever should be done, like replicating
      ## old phenoData variables length(reporters) times
      warning("Unexpected number of samples in phenoData slot. Dropping it.")
    }
  }
  phenoData(msnset) <- .phenoData
  ## Updating protocol slot 
  if (nrow(protocolData(object))>0) { 
    if (nrow(protocolData(object))==length(reporters)) {
      protocolData(msnset) <- protocolData(object)
    } else {
      ## Here, something more clever should be done, like replicating
      ## old phenoData variables length(reporters) times
      warning("Unexpected number of features in featureData slot. Dropping it.")
    }
  }
  ## Returning shiny MSnSet object
  if (validObject(msnset))
    return(msnset)
}


