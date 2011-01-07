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
                          tic=unlist(eapply(spectra(object),
                            function(x) sum(intensity(x)))),
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
  object@process@processing <- c(object@process@processing,
                                 paste(n,"precursors extracted:",date()))
  ## TODO: need to update phenoData - check files
  return(new("MSnExp",
             spectra=list2env(object@spectra[sel]),
             process=object@process,
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
  object@process@removedPeaks <- c(object@process@removedPeaks,
                                   as.character(t))
  object@process@processing <- c(object@process@processing,
                                 paste("Curves <= ",t," set to '0': ",date(),sep=""))
  return(object)
}


bg.correct.MSnExp <- function(object,bg,verbose=TRUE) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <-  llply(object@spectra,function(x) bg.correct(x,bg),.progress=progress)
  object@spectra <- spectraList
  object@process@removedPeaks <- c(object@process@removedPeaks,
                                   as.character(bg))
  if (bg<0)
    bg <- "min. int."
  object@process@processing <- c(object@process@processing,
                                 paste("Backgroung corrected using ",bg,": ",date(),sep=""))
  return(object)
}


clean.MSnExp <- function(object,verbose=TRUE) {
  pb <- txtProgressBar(min = 0, max = total, style = 3) 
  ifelse(verbose,progress <- "text",progress <- "none")
  spectra <- llply(spectra(object),function(x) clean(x),.progress=progress)
  object@assayData <- list2env(spectra)
  object@process@cleaned <- TRUE
  object@process@processing <- c(object@process@processing,
                                 paste("Spectra cleaned: ",date(),sep=""))
  return(object)
}

quantify.MSnExp <- function(object,reporters,method,verbose) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <- spectra(object)
  ## Quantification -- creating exprs for assayData slot
  if (length(spectraList)==1) {
    peakData <- quantify(spectraList[[1]],reporters,method)
    .exprs <- t(peakData$peakArea)
    .qual <- t(peakData$curveData)
  } else {
    peakData <- laply(spectraList,quantify,reporters,method,.progress=progress)
    .exprs <- do.call(rbind,sapply(peakData,"[","peakArea"))
    .qual <- do.call(rbind,sapply(peakData,"[","curveStats"))
  }
  prec <- sapply(spectraList,precursorMz)
  feat <- make.unique(as.character(prec))
  rownames(.qual) <- 1:nrow(.qual)
  rownames(.exprs) <- feat
  colnames(.exprs) <- reporters@reporterNames
  ## Updating MSnprocess slot
  object@process@processing <- c(object@process@processing,
                                 paste(reporters@name," quantification by ",method,
                                       ": ",date(),sep=""))
  object@process@centroided <- TRUE
  ## Creating new MSnSet
  msnset <- new("MSnSet",
                qual=.qual,
                exprs=.exprs, 
                process=object@process,
                files=object@files,
                protocolData=protocolData(object),
                experimentData=experimentData(object))
  ## Updating featureData slot or creating one
  fd <- header(object)
  if (nrow(fData(object))>0) { 
    if (nrow(fData(object))==length(object)) {
      .featureData <- new("AnnotatedDataFrame",data=cbind(fData(object),fd))
    } else {
      warning("Unexpected number of features in featureData slot. Dropping it.")
    }
  } else {
    rownames(fd) <- feat
    .featureData <- new("AnnotatedDataFrame",data=fd)
  }
  featureData(msnset) <- .featureData
  ## Updating phenoData slot or creating one
  pd <- data.frame(mz=reporters@mz,
                   reporters=reporters@name,
                   row.names=reporters@reporterNames)
  if (nrow(pData(object))>0) { 
    if (nrow(pData(object))==length(reporters)) {
      .phenoDataData <- new("AnnotatedDataFrame",data=cbind(pData(object),pd))
    } else {
      warning("Unexpected number of samples in phenoData slot. Dropping it.")
    }
  } else {
    .phenoData <- new("AnnotatedDataFrame",data=pd)
  }
  phenoData(msnset) <- .phenoData
  ## Updating protocol slot
  if (nrow(protocolData(object))>0) { 
    if (nrow(protocolData(object))==length(reporters)) 
      protocolData(msnset) <- protocolData(object)
    else 
      warning("Unexpected number of features in featureData slot. Dropping it.")
  }
  ## Returning shiny MSnSet object
  if (validObject(msnset))
    return(msnset)
}


