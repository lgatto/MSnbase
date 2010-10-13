show.MSnExp <- function(object) {
  cat("Object of class \"",class(object),"\"\n",sep="")
  cat(" Object size in memory: ")
  print(object.size(object),units="Mb")
  cat("- - - Meta data  - - -\n")
  cat(" Data description:",object@description,"\n")  
  cat(" Loaded from:\n")
  for (i in 1:length(object@files)) {
    f <- basename(object@files[i])
    cat("   ",f,"\n")
  }
  show(object@process)
  cat("- - - Spectra data - - -\n")
  msnLevels <- unique(msLevel(object))
  cat(" MSn level(s):",msnLevels,"\n")
  if (all(msLevel(object)>1)) {
    cat(" Number of MS1 acquisitions:",length(unique(ms1scanNum(object))),"\n")
    cat(" Number of MS2 scans:",length(spectra(object)),"\n")
    msnPrecMz <- precursorMz(object)
    nbPrecIons <- length(msnPrecMz)
    cat(" Number of precursor ions:",nbPrecIons,"\n")
    if (nbPrecIons>0) {
      cat("",length(unique(msnPrecMz)),"unique MZs\n")
      cat(" Precursor MZ's:",paste(signif(range(msnPrecMz),5),collapse=" - "),"\n")
    }
    msnMzRange <- round(range(mz(object)),2)
    cat(" MSn M/Z range:",msnMzRange,"\n")
  } else {
    cat(" Number of MS1 scans:",length(spectra(object)),"\n")
  }
  msnRt <- rtime(object)
  if (length(msnRt)>0) {
    rtr <- range(msnRt)
    cat(" MSn retention times:",formatRt(rtr[1]),"-",formatRt(rtr[2]),"minutes\n")
  }
}

"[.MSnExp" <- function(x,i) {
  if (max(i)>length(x) | min(i)<1)
    stop("subscript out of bonds")
  if (length(i)==1)
    return(spectra(x)[[i]])
  return(spectra(x)[i])
}

header.MSnExp <- function(object) {
  if (any(msLevel(object)<2))
    stop("header() only works for MS levels > 1.")
  tbl <- table(object@fromFile)
  idx <- as.numeric(unlist(apply(tbl,1,function(x) 1:x)))
  return(data.frame(cbind(index=idx,
                          file=object@fromFile,
                          retention.time=rtime(object),
                          precursor.mz=precursorMz(object),
                          peaks.count=peaksCount(object),
                          tic=sapply(spectra(object),
                            function(x) sum(intensity(x))),
                          ms.level=msLevel(object),
                          charge=precursorCharge(object),
                          collision.energy=collisionEnergy(object))))
       }

mergeSpectra <- function(object, ## MSnExp object
                         fun=sum,
                         verbose=TRUE) {
  spectra <- object@spectra
  prec <- sapply(spectra,function(x) x@precursorMz)
  uprec <- unique(prec)
  luprec <- length(uprec)
  if (verbose) 
    cat(luprec,"unique precursors for",length(spectra),"spectra\n")
  newSpectra <- vector("list",length=luprec)
  if (verbose)
    pb <- txtProgressBar(min = 0, max = luprec, style = 3)
  for (i in 1:luprec) {
    if (verbose)
      setTxtProgressBar(pb, i)   
    pi <- uprec[i]
    sel <- prec %in% pi
    l <- spectra[sel]
    ## unique M/Z values
    mzs <- unique(unlist(lapply(l,function(x) x@mz)))
    ## vectors of intensitues (ints) and number of times
    ## a given itensity is found (icounts)
    ints <- icounts <- numeric(length(mzs))
    names(icounts) <- mzs
    allints <- unlist(lapply(l,function(x) x@intensity))
    names(allints) <- unlist(lapply(l,function(x) x@mz))
    for (j in 1:length(mzs)) {
      k <- as.character(mzs[j])
      ints[j] <- fun(allints[names(allints) %in% k])
      icounts[j] <- sum(names(allints) %in% k)
    }
    o <- order(mzs)
    newSpectra[[i]] <- new("Spectrum",
                           merged=which(sel),
                           rt=unique(sapply(l,function(x) x@rt)),
                           msLevel=unique(sapply(l,function(x) x@msLevel)),
                           precursorMz=pi,
                           peaksCount=length(ints),
                           intensity=ints[o],
                           mz=mzs[o])
    validObject(newSpectra[[i]])
  }
  if (verbose)
    close(pb)
  object@process@processing <- c(object@process@processing,
                                 paste("Precursor ions with identical M/Z merged:",
                                       date()))
  object@process@merged <- TRUE
  return(new("MSnExp",
             spectra=newSpectra,
             metadata=object@metadata,
             process=object@process))
}

extractPrecSpectra <- function(object,prec) {
  sel <- sapply(object@spectra,function(x) x@precursorMz %in% prec)
  n <- length(prec)
  object@process@processing <- c(object@process@processing,
                                 paste(n,"precursors extracted:",date()))
  return(new("MSnExp",
             spectra=object@spectra[sel],
             process=object@process,
             description=object@description,
             fromFile=object@fromFile[sel],
             files=object@files,
             assayData=object@assayData,
             phenoData=object@phenoData,
             featureData=object@featureData,
             experimentData=object@experimentData,
             annotation=object@annotation,
             protocolData=object@protocolData))
}


extractSpectra <- function(object,selected) {
  if (!is.logical(selected))
    stop("Please use logicals to extract a set of spectra.")
  object@spectra <- spectra(object)[selected]
  object@fromFile <- object@fromFile[selected]
  return(object)
}

removePeaks.MSnExp <- function(object,t,verbose=TRUE) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <-  llply(object@spectra,function(x) removePeaks(x,t),.progress=progress)
  object@spectra <- spectraList
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
  ifelse(verbose,progress <- "text",progress <- "none")
  spectra <- llply(object@spectra,function(x) clean(x),.progress=progress)
  object@spectra <- spectra
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
                                 paste("Quantification by ",method,
                                       reporters@name,": ",date(),sep=""))
  object@process@centroided <- TRUE
  ## Creating new MSnSet
  msnset <- new("MSnSet",
                qual=.qual,
                exprs=.exprs, 
                process=object@process,
                description=object@description,
                files=object@files,
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
  return(msnset)
}


