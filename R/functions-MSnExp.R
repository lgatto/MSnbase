## mergeSpectra <- function(object, ## MSnExp object
##                          fun = sum,
##                          verbose = isMSnbaseVerbose()) {
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

removePeaks_MSnExp <- function(object, t = "min", verbose = isMSnbaseVerbose()) {
  ifelse(verbose, progress <- "text", progress <- "none")
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


clean_MSnExp <- function(object, all, verbose = isMSnbaseVerbose()) {
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
    lockEnvironment(e, bindings = TRUE)
    object@assayData <- e
    if (validObject(object))
        return(object)
}


normalise_MSnExp <- function(object,method) {
    ## Can not directly assign to assayData, as that environment is locked!
    e <- new.env()
    sapply(featureNames(object),
           function(x) {
               sp <- get(x,envir=assayData(object))
               xx <- normalise(sp,method)
               assign(x,xx,envir=e)
               invisible(TRUE)
           })
    object@processingData@processing <- c(object@processingData@processing,
                                          paste0("Spectra normalised (",method,"): ",
                                                 date()))
    object@processingData@normalised <- TRUE
    lockEnvironment(e, bindings = TRUE)
    object@assayData <- e
    if (validObject(object))
        return(object)
}

bin_MSnExp <- function(object, binSize = 1, verbose = isMSnbaseVerbose()) {
    ## copied from clean_MSnExp
    e <- new.env()

    if (verbose) {
        ._cnt <- 1
        pb <- txtProgressBar(min = 0, max = length(object), style = 3)
    }

    mzrange <- range(eapply(assayData(object), mz))
    breaks <- seq(floor(mzrange[1]), ceiling(mzrange[2]), by = binSize)

    sapply(featureNames(object),
           function(x) {
               if (verbose) {
                   setTxtProgressBar(pb, ._cnt)
                   ._cnt <<- ._cnt+1
               }
               sp <- get(x, envir = assayData(object))
               xx <- bin_Spectrum(sp, breaks = breaks)
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
    lockEnvironment(e, bindings = TRUE)
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
                             ..., verbose = isMSnbaseVerbose()) {
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
    lockEnvironment(e, bindings = TRUE)    
    object@assayData <- e
    if (validObject(object))
    return(object)
}

smooth_MSnExp <- function(object, method, halfWindowSize, ...,
                          verbose = isMSnbaseVerbose()) {
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
    lockEnvironment(e, bindings = TRUE)    
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

removeReporters_MSnExp <- function(object, reporters = NULL,
                                   clean = FALSE, verbose = isMSnbaseVerbose()) {
  ifelse(verbose, progress <- "text", progress <- "none")
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

#' @title Estimate m/z scattering in consecutive scans
#'
#' @description
#'
#' Estimate scattering of m/z values (due to technical, instrument specific
#' noise) for the same ion in consecutive scans of a LCMS experiment.
#'
#' @details
#'
#' The m/z values of the same ions in consecutive scans (spectra) of a LCMS run
#' will not be identical. This random noise is expected to be smaller than the
#' resolution of the MS instrument. The distribution of differences of m/z
#' values from neighboring spectra is thus expected to be (at least) bi-modal
#' with the first peak representing the above described random variation and
#' the second (or largest) peak the m/z resolution. The m/z value of the first
#' local minimum between these first two peaks in the distribution is returned
#' as the *m/z scattering*.
#'
#' @param x `MSnExp` or `OnDiskMSnExp` object.
#'
#' @param halfWindowSize `integer(1)` defining the half window size for the
#'     moving window to combine consecutive spectra.
#' 
#' @author Johannes Rainer
#'
#' @md
#'
#' @seealso [estimateMzResolution()] for the function to estimate a
#'     profile-mode spectrum's m/z resolution from it's data.
#'
#' @examples
#' 
#' library(MSnbase)
#' library(msdata)
#' ## Load the MS1 data of a profile-mode file
#' f <- msdata::proteomics(full.names = TRUE,
#'     pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
#'
#' im <- readMSData(f, mode = "inMem", msLevel = 1L)
#'
#' res <- estimateMzScattering(im)
#'
#' ## Plot the distribution of estimated m/z scattering
#' plot(density(unlist(res)))
#'
#' ## Compare the m/z resolution and m/z scattering of the spectrum with the
#' ## most peaks
#' idx <- which.max(unlist(spectrapply(im, peaksCount)))
#'
#' res[[idx]]
#' estimateMzResolution(im[[idx]])
#' ## As expected, the m/z scattering is much lower than the m/z resolution.
estimateMzScattering <- function(x, halfWindowSize = 1L) {
    if (!is(x, "MSnExp"))
        stop("'x' should be a 'MSnExp' object")
    res <- lapply(split(spectra(x), fromFile(x)), function(z) {
        len_z <- length(z)
        mzs <- vector("list", len_z)
        for (i in seq_along(z)) {
            mzs[[i]] <- .estimate_mz_scattering(
                sort(unlist(lapply(z[max(1, i - halfWindowSize):
                                     min(i + halfWindowSize, len_z)], mz))))
        }
        mzs
    })
    res <- unlist(res, fromFile(x))
    names(res) <- featureNames(x)
    res
}

#' @title Combine signal from consecutive spectra of LCMS experiments
#'
#' @description
#'
#' `combineSpectraMovingWindow` combines signal from consecutive spectra within
#' a file. The resulting `MSnExp` has the same total number of spectra than the
#' original object, but with each individual's spectrum information
#' representing aggregated data from the original spectrum and its neighboring
#' spectra.
#'
#' Note that the function returns always a `MSnExp` object, even if `x` was an
#' `OnDiskMSnExp` object.
#' 
#' @details
#'
#' The method assumes same ions being measured in consecutive scans (i.e. LCMS
#' data) and thus combines their signal which can increase the increase the
#' signal to noise ratio.
#' 
#' See [combineSpectra()] for details.
#' 
#' @param x `MSnExp` or `OnDiskMSnExp` object.
#'
#' @param halfWindowSize `integer(1)` with the half window size for the moving
#'     window.
#' 
#' @param BPPARAM parallel processing settings.
#' 
#' @inheritParams combineSpectra
#'
#' @return `MSnExp` with the same number of spectra than `x`.
#'
#' @md
#'
#' @seealso
#'
#' [combineSpectra()] for the function combining spectra provided in
#' a `list`.
#'
#' [estimateMzScattering()] for a function to estimate m/z value scattering in
#' consecutive spectra.
#' 
#' @author Johannes Rainer
combineSpectraMovingWindow <- function(x, halfWindowSize = 1L,
                                       mzFun = base::mean,
                                       intensityFun = base::sum, mzd,
                                       BPPARAM = bpparam()){
    if (!is(x, "MSnExp"))
        stop("'x' has to be a 'MSnExp' or an 'OnDiskMSnExp'")
    if (is(x, "OnDiskMSnExp"))
        x <- as(x, "MSnExp")
    ## Combine spectra per file
    new_sp <- bplapply(split(spectra(x), fromFile(x)), function(z, intF,
                                                                mzF, hws, mzd) {
        len_z <- length(z)
        res <- vector("list", len_z)
        for (i in seq_along(z)) {
            res[[i]] <- combineSpectra(z[max(1, i - hws):min(i + hws, len_z)],
                                       mzFun = mzF, intensityFun = intF,
                                       main = hws + 1L, mzd = mzd)
        }
        res
    }, intF = intensityFun, mzF = mzFun, hws = as.integer(halfWindowSize),
    mzd = mzd, BPPARAM = BPPARAM)
    new_sp <- unsplit(new_sp, fromFile(x))
    names(new_sp) <- featureNames(x)
    x@assayData <- list2env(new_sp)
    lockEnvironment(x@assayData, bindings = TRUE)
    if (validObject(x))
        x
}
