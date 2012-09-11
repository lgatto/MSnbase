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

extractPrecSpectra.MSnExp <- function(object,prec) {
  nmissing <- sum(!prec %in% precursorMz(object))
  if (nmissing!=0)
    warning(nmissing," precursor MZ values not found in 'MSnExp' object.") 
  sel <-precursorMz(object) %in% prec
  orghd <- header(object)
  nms <- names(precursorMz(object)[sel])
  n <- length(prec)
  m <- length(nms)
  ## updating object
  object@processingData@processing <- c(object@processingData@processing,
                                        paste(n,"(",m,
                                              ") precursors (spectra) extracted: ",
                                              date(),sep=""))
  object@assayData <- list2env(mget(nms,assayData(object)))
  object@featureData <- object@featureData[nms,]
  if (object@.cache$level > 0)
    object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                      hd = orghd[sel, ]),
                                 object@.cache$level)
  if (validObject(object))
    return(object)
}


## Defunct in v 1.5.2
## extractSpectra.MSnExp <- function(object,selected) {
##   slist <- spectra(object)[selected]
##   object@assayData <- list2env(slist)
##   object@featureData <- object@featureData[selected,]
##   object@processingData@processing <- c(object@processingData@processing,
##                                         paste(sum(selected),
##                                               " spectra extracted: ",
##                                               date(),sep=""))
##   if (object@.cache$level > 0)
##     object@.cache <- setCacheEnv(assayData(object), object@.cache$level)
##   if (validObject(object))
##     return(object)
## }

removePeaks.MSnExp <- function(object,t="min",verbose=TRUE) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <-  llply(spectra(object),
                        function(x) removePeaks(x,t),
                        .progress = progress)
  object@assayData <- list2env(spectraList)
  object@processingData@removedPeaks <- c(object@processingData@removedPeaks,
                                          as.character(t))
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Curves <= ",
                                              t
                                              ," set to '0': ",
                                              date(),sep=""))
  if (object@.cache$level > 0) {
    hd <- header(object)
    hd$ionCount <- ionCount(object)
    object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                      hd = hd),
                                 object@.cache$level)
  }
  if (validObject(object))
    return(object)
}


clean.MSnExp <- function(object,verbose=TRUE) {
  ## -- was ---------------------------------------------------
  ##  ifelse(verbose,progress <- "text",progress <- "none")
  ##  spectra <- llply(spectra(object),function(x) clean(x),.progress=progress)
  ##  object@assayData <- list2env(spectra)
  ## -- new ---------------------------------------------------
  if (verbose) {
    ._cnt <- 1
    pb <- txtProgressBar(min = 0, max = length(object), style = 3)
  }
  sapply(featureNames(object),
         function(x) {
           if (verbose) {
             setTxtProgressBar(pb, ._cnt)
             ._cnt <<- ._cnt+1
           }
           sp <- get(x,envir=assayData(object))
           xx <- clean(sp)
           assign(x,xx,envir=assayData(object))
           invisible(TRUE)
         })
  if (verbose) {
    close(pb)
    rm(pb)
    rm(._cnt)
  }
  ## ----------------------------------------------------------
  object@processingData@cleaned <- TRUE
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Spectra cleaned: ",date(),sep=""))
  
  if (object@.cache$level > 0) {
    hd <- header(object)
    hd$peaks.count <- peaksCount(object)
    object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                      hd = hd),
                                 object@.cache$level)
  }
  if (validObject(object))
    return(object)
}

quantify.MSnExp <- function(object, method, reporters, strict, parallel, verbose) {
  ## Display progress bar with eapply
  ## TODO - test if using eapply is more efficient in terms of mem/cpu usage
  ## if (verbose) {
  ##   ._cnt <- 1
  ##   pb <- txtProgressBar(min = 0, max = length(object), style = 3)
  ##   ## Quantification -- creating exprs for assayData slot
  ##   peakData <- eapply(assayData(object),function(x) {
  ##     setTxtProgressBar(pb, ._cnt)
  ##     ._cnt <<- ._cnt+1
  ##     quantify(x,method,reporters)
  ##   })
  ##   close(pb)
  ##   rm(pb)
  ##   rm(._cnt)
  ## } else {
  ##   peakData <- eapply(assayData(object),quantify,method,reporters)
  ## }
  ifelse(verbose, progress <- "text", progress <- "none")
  if (.Platform$OS.type == "windows") {
    parallel <- FALSE
    if (verbose)
      message("Parallel processing not yet supported on Windows.")    
  }
  spectraList <- spectra(object)
  ## Quantification -- creating exprs for assayData slot
  if (length(spectraList) == 1) {
    peakData <- quantify(spectraList[[1]], method, reporters, strict)
    .exprs <- t(peakData$peakQuant)
    .qual <- t(peakData$curveData)
  } else {
    if (parallel & require(foreach) & require(doMC)) {
      registerDoMC()
    }
    peakData <- llply(spectraList, quantify, method, reporters, strict,
                      .progress = progress, .parallel = parallel)      
    .exprs <- do.call(rbind, sapply(peakData, "[", "peakQuant"))
    ## .qual <- do.call(rbind, sapply(peakData, "[", "curveStats")) ## Time consuming - consider removing or caching
    .qual <- data.frame()
  }
  rownames(.exprs) <- sub(".peakQuant", "", rownames(.exprs))  
  rownames(.qual) <- sub(".curveStats", "", rownames(.qual)) 
  ## Updating MSnprocess slot
  object@processingData@processing <- c(object@processingData@processing,
                                        paste(reporters@name,
                                              ifelse(strict, " (strict) ", " "),
                                              "quantification by ", method,
                                              ": ", date(), sep=""))
  ## Creating new featureData slot or creating one
  if (verbose)
    message("Preparing meta-data")
  fd <- header(object) ## Time consuming - consider caching
  if (nrow(fData(object)) > 0) { 
    if (nrow(fData(object)) == length(object)) {
      fd <- combine(fData(object), fd)
    } else {
      warning("Unexpected number of features in featureData slot. Dropping it.")
    }
  }
  ## featureData rows must be reordered to match assayData rows
  .featureData <- new("AnnotatedDataFrame", data=fd[rownames(.exprs), ])
  ## Creating new phenoData slot or creating one
  .phenoData <- new("AnnotatedDataFrame",
                    data = data.frame(mz = reporters@mz,
                      reporters = reporters@name,
                      row.names = reporters@reporterNames))
  if (nrow(pData(object)) > 0) { 
    if (nrow(pData(object)) == length(reporters)) {
      .phenoData <- combine(phenoData(object), .phenoData)
    } else {
      ## Here, something more clever should be done, like replicating
      ## old phenoData variables length(reporters) times
      if (verbose)
        message("Original MSnExp and new MSnSet have different number of samples in phenoData. Dropping original.")
    }
  }
  if (verbose)
    message("Creating 'MSnSet' object")
  msnset <- new("MSnSet",
                qual = .qual,
                exprs = .exprs, 
                experimentData = experimentData(object),
                phenoData = .phenoData,
                featureData = .featureData,
                annotation = "No annotation")
  
  ## copying processingData  
  msnset@processingData <- object@processingData
  
  ## Updating protocol slot 
  if (nrow(protocolData(object)) > 0) { 
    if (nrow(protocolData(object)) == length(reporters)) {
      .protocolData <- protocolData(object)
    } else {
      warning("protocolData does not match with reporters. Dropping it.")
    }
  }
  ## Returning shiny MSnSet object
  if (validObject(msnset))
    return(msnset)
}


normalise.MSnExp <- function(object,method) {
  sapply(featureNames(object),
         function(x) {
           sp <- get(x,envir=assayData(object))
           xx <- normalise(sp,method)
           assign(x,xx,envir=assayData(object))
           invisible(TRUE)
         })
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Spectra normalised (",method,"): ",
                                              date(),
                                              sep=""))
  object@processingData@normalised <- TRUE
  if (validObject(object))
    return(object)
}

precSelection <- function(object,n=NULL) {
  allPrecs <- precursorMz(object)
  if (!is.null(n))
    allPrecs <- round(allPrecs,n)
  number.selection <- c()
  scanNums <- MSnbase:::precScanNum(object)
  uprecmz <- unique(allPrecs)
  for (mp in uprecmz)
      number.selection <- c(number.selection,
                            length(unique(scanNums[allPrecs==mp])))
  names(number.selection) <- uprecmz
  return(number.selection)
}

precSelectionTable <- function(object,...) {
  x <- precSelection(object,...)
  return(table(x))
}

removeReporters.MSnExp <- function(object,reporters=NULL,clean=FALSE,verbose=TRUE) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <-  llply(spectra(object),function(x) removeReporters(x,reporters,clean),.progress=progress)
  object@assayData <- list2env(spectraList)
  repname <- names(reporters)
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Removed", repname, "reporter ions",sep=" "))
  if (validObject(object))
    return(object)
}
