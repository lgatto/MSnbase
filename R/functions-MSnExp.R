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

extractPrecSpectra_MSnExp <- function(object,prec) {
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
                                        paste(n," (",m,
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

removePeaks_MSnExp <- function(object,t="min",verbose=TRUE) {
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


clean_MSnExp <- function(object, all, verbose = TRUE) {
  ## -- was ---------------------------------------------------
  ##  ifelse(verbose,progress <- "text",progress <- "none")
  ##  spectra <- llply(spectra(object),function(x) clean(x),.progress=progress)
  ##  object@assayData <- list2env(spectra)
  ## -- new ---------------------------------------------------
  e <- new.env()
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
           sp <- get(x, envir = assayData(object))
           xx <- clean(sp, all)
           assign(x, xx, envir = e)
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
                                        paste0("Spectra cleaned: ", date()))

  if (object@.cache$level > 0) {
    hd <- header(object)
    hd$peaks.count <- peaksCount(object)
    object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                      hd = hd),
                                 object@.cache$level)
  }
  object@assayData <- e
  if (validObject(object))
    return(object)
}


normalise_MSnExp <- function(object,method) {
  sapply(featureNames(object),
         function(x) {
           sp <- get(x,envir=assayData(object))
           xx <- normalise(sp,method)
           assign(x,xx,envir=assayData(object))
           invisible(TRUE)
         })
  object@processingData@processing <- c(object@processingData@processing,
                                        paste0("Spectra normalised (",method,"): ",
                                              date()))
  object@processingData@normalised <- TRUE
  if (validObject(object))
    return(object)
}

bin_MSnExp <- function(object, binSize=1, verbose=TRUE) {
  ## copied from clean_MSnExp
  e <- new.env()

  if (verbose) {
    ._cnt <- 1
    pb <- txtProgressBar(min = 0, max = length(object), style = 3)
  }

  mzrange <- range(eapply(assayData(object), mz))
  breaks <- seq(floor(mzrange[1]), ceiling(mzrange[2]), by=binSize)

  sapply(featureNames(object),
         function(x) {
           if (verbose) {
             setTxtProgressBar(pb, ._cnt)
             ._cnt <<- ._cnt+1
           }
           sp <- get(x, envir = assayData(object))
           xx <- bin(sp, breaks=breaks)
           assign(x, xx, envir = e)
           invisible(TRUE)
         })
  if (verbose) {
    close(pb)
    rm(pb)
    rm(._cnt)
  }
  ## ----------------------------------------------------------
  object@processingData@processing <- c(object@processingData@processing,
                                        paste0("Spectra binned: ", date()))
  if (object@.cache$level > 0) {
    hd <- header(object)
    hd$peaks.count <- peaksCount(object)
    object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                      hd = hd),
                                 object@.cache$level)
  }
  object@assayData <- e
  if (validObject(object))
    return(object)
}

compare_MSnExp <- function(object, fun, ...) {

  nm <- featureNames(object)
  cb <- combn(nm, 2, function(x) {
    compare_Spectra(object[[x[1]]], object[[x[2]]], fun=fun, ...)
  })
  m <- matrix(NA, length(object), length(object),
              dimnames=list(nm, nm))
  ## fill lower triangle of the matrix
  m[lower.tri(m)] <- cb
  ## copy to upper triangle
  for (i in 1:nrow(m)) {
    m[i, ] <- m[, i]
  }

  return(m)
}

pickPeaks_MSnExp <- function(object, halfWindowSize, method, SNR,
                             ..., verbose = TRUE) {
  ## copied from clean_MSnExp
  e <- new.env()

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
           sp <- get(x, envir = assayData(object))
           xx <- pickPeaks(sp, halfWindowSize = halfWindowSize,
                           method = method, SNR = SNR, ...)
           assign(x, xx, envir = e)
           invisible(TRUE)
         })
  if (verbose) {
    close(pb)
    rm(pb)
    rm(._cnt)
  }
  ## ----------------------------------------------------------
  ## TODO: @lgatto object@processingData@centroided ?
  ## object@processingData@smoothed <- TRUE
  object@processingData@processing <- c(object@processingData@processing,
                                        paste0("Spectra centroided: ",
                                               date()))
  if (object@.cache$level > 0) {
    hd <- header(object)
    hd$peaks.count <- peaksCount(object)
    object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                      hd = hd),
                                 object@.cache$level)
  }
  object@assayData <- e
  if (validObject(object))
    return(object)
}

smooth_MSnExp <- function(object, method, halfWindowSize, ..., verbose = TRUE) {
  ## copied from clean_MSnExp
  e <- new.env()

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
           sp <- get(x, envir = assayData(object))
           xx <- smooth(sp, method = method, halfWindowSize = halfWindowSize,
                        ...)
           assign(x, xx, envir = e)
           invisible(TRUE)
         })
  if (verbose) {
    close(pb)
    rm(pb)
    rm(._cnt)
  }
  ## ----------------------------------------------------------
  object@processingData@smoothed <- TRUE
  object@processingData@processing <- c(object@processingData@processing,
                                        paste0("Spectra smoothed (",
                                               method, "): ", date()))
  if (object@.cache$level > 0) {
    hd <- header(object)
    hd$peaks.count <- peaksCount(object)
    object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                      hd = hd),
                                 object@.cache$level)
  }
  object@assayData <- e
  if (validObject(object))
    return(object)
}


precSelection <- function(object,n=NULL) {
  allPrecs <- precursorMz(object)
  if (!is.null(n))
    allPrecs <- round(allPrecs,n)
  number.selection <- c()
  scanNums <- precScanNum(object)
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

removeReporters_MSnExp <- function(object, reporters=NULL, clean=FALSE, verbose=TRUE) {
  ifelse(verbose,progress <- "text",progress <- "none")
  spectraList <-  llply(spectra(object),
                        function(x) removeReporters(x, reporters, clean),
                        .progress = progress)
  object@assayData <- list2env(spectraList)
  repname <- names(reporters)
  object@processingData@processing <- c(object@processingData@processing,
                                        paste("Removed", repname, "reporter ions",sep=" "))
  if (validObject(object))
    return(object)
}

