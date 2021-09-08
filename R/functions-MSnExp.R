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
                             ..., msLevel., verbose = isMSnbaseVerbose()) {
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
                           method = method, SNR = SNR, msLevel. = msLevel.,
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

smooth_MSnExp <- function(object, method, halfWindowSize, ..., msLevel.,
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
                        msLevel. = msLevel., ...)
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
#' @note
#'
#' For `timeDomain = TRUE` the function does **not** return the estimated
#' scattering of m/z values, but the scattering of `sqrt(mz)` values.
#'
#' @param x `MSnExp` or `OnDiskMSnExp` object.
#'
#' @param halfWindowSize `integer(1)` defining the half window size for the
#'     moving window to combine consecutive spectra.
#'
#' @param timeDomain `logical(1)` whether m/z scattering should be estimated
#'     on `mz` (`timeDomain = FALSE`) or `sqrt(mz)` (`timeDomain = TRUE`)
#'     values. See [combineSpectraMovingWindow()] for details on this
#'     parameter.
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
#' ## Load a profile-mode LC-MS data file
#' f <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)[1]
#' od <- readMSData(f, mode = "onDisk")
#' im <- as(filterRt(od, c(10, 20)), "MSnExp")
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
#' abline(v = res[[idx]], lty = 2)
#' estimateMzResolution(im[[idx]])
#' ## As expected, the m/z scattering is much lower than the m/z resolution.
estimateMzScattering <- function(x, halfWindowSize = 1L, timeDomain = FALSE) {
    if (!is(x, "MSnExp"))
        stop("'x' should be a 'MSnExp' object")
    res <- lapply(split(spectra(x), fromFile(x)),
                  FUN = .estimate_mz_scattering_list,
                  halfWindowSize = halfWindowSize, timeDomain = timeDomain)
    res <- unsplit(res, fromFile(x))
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
#' spectra. This is thus equivalent with a smoothing of the data in retention
#' time dimension.
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
#' Intensities (and m/z values) for signals with the same m/z value in
#' consecutive scans are aggregated using the `intensityFun`.
#' m/z values of intensities from consecutive scans will never be exactly
#' identical, even if they represent signal from the same ion. The function
#' determines thus internally a similarity threshold based on differences
#' between m/z values within and between spectra below which m/z values are
#' considered to derive from the same ion. For robustness reasons, this
#' threshold is estimated on the 100 spectra with the largest number of
#' m/z - intensity pairs (i.e. mass peaks).
#'
#' See [meanMzInts()] for details.
#'
#' Parameter `timeDomain`: by default, m/z-intensity pairs from consecutive
#' scans to be aggregated are defined based on the square root of the m/z
#' values. This is because it is highly likely that in all QTOF MS instruments
#' data is collected based on a timing circuit (with a certain variance) and
#' m/z values are later derived based on the relationship `t = k * sqrt(m/z)`.
#' Differences between individual m/z values will thus be dependent on the
#' actual m/z value causing both the difference between m/z values and their
#' scattering being different in the lower and upper m/z range. Determining
#' m/z values to be combined on the `sqrt(mz)` reduces this dependency. For
#' non-QTOF MS data `timeDomain = FALSE` might be used instead.
#'
#' @note
#'
#' The function has to read all data into memory for the spectra combining
#' and thus the memory requirements of this function are high, possibly
#' preventing its usage on large experimental data. In these cases it is
#' suggested to perform the combination on a per-file basis and save the
#' results using the [writeMSData()] function afterwards.
#'
#' @param x `MSnExp` or `OnDiskMSnExp` object.
#'
#' @param halfWindowSize `integer(1)` with the half window size for the moving
#'     window.
#'
#' @param ppm `numeric(1)` to define an m/z relative deviation. Note that if
#'     only `ppm` should be considered but not `mzd`, `mzd` should be set to
#'     `0` (i.e. `mzd = 0`). This parameter is directly passed to
#'     [meanMzInts()].
#'
#' @param BPPARAM parallel processing settings.
#'
#' @inheritParams meanMzInts
#'
#' @return `MSnExp` with the same number of spectra than `x`.
#'
#' @md
#'
#' @seealso
#'
#' [meanMzInts()] for the function combining spectra provided in
#' a `list`.
#'
#' [estimateMzScattering()] for a function to estimate m/z value scattering in
#' consecutive spectra.
#'
#' @author Johannes Rainer, Sigurdur Smarason
#'
#' @examples
#'
#' library(MSnbase)
#' library(msdata)
#'
#' ## Read a profile-mode LC-MS data file.
#' fl <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)[1]
#' od <- readMSData(fl, mode = "onDisk")
#'
#' ## Subset the object to the retention time range that includes the signal
#' ## for proline. This is done for performance reasons.
#' rtr <- c(165, 175)
#' od <- filterRt(od, rtr)
#'
#' ## Combine signal from neighboring spectra.
#' od_comb <- combineSpectraMovingWindow(od)
#'
#' ## The combined spectra have the same number of spectra, same number of
#' ## mass peaks per spectra, but the signal is larger in the combined object.
#' length(od)
#' length(od_comb)
#'
#' peaksCount(od)
#' peaksCount(od_comb)
#'
#' ## Comparing the chromatographic signal for proline (m/z ~ 116.0706)
#' ## before and after spectra data combination.
#' mzr <- c(116.065, 116.075)
#' chr <- chromatogram(od, rt = rtr, mz = mzr)
#' chr_comb <- chromatogram(od_comb, rt = rtr, mz = mzr)
#'
#' par(mfrow = c(1, 2))
#' plot(chr)
#' plot(chr_comb)
#' ## Chromatographic data is "smoother" after combining.
combineSpectraMovingWindow <- function(x, halfWindowSize = 1L,
                                       intensityFun = base::mean,
                                       mzd = NULL,
                                       timeDomain = FALSE,
                                       weighted = FALSE,
                                       ppm = 0,
                                       BPPARAM = bpparam()){
    if (!is(x, "MSnExp"))
        stop("'x' has to be a 'MSnExp' or an 'OnDiskMSnExp'")
    if (is(x, "OnDiskMSnExp"))
        x <- as(x, "MSnExp")
    ## Combine spectra per file
    new_sp <- bplapply(split(spectra(x), fromFile(x)), FUN = function(z, intF,
                                                                      wght, hws,
                                                                      mzd,
                                                                      timeD,
                                                                      ppm) {
        len_z <- length(z)
        ## Estimate m/z scattering on the 100 spectra with largest number of
        ## peaks
        if (is.null(mzd)) {
            idx <- order(unlist(lapply(z, function(y) y@peaksCount)),
                         decreasing = TRUE)[1:min(100, len_z)]
            mzs <- .estimate_mz_scattering_list(z[idx], halfWindowSize = hws,
                                                timeDomain = timeD)
            dens <- .density(unlist(mzs))
            mzd <- dens$x[which.max(dens$y)]
        }
        ## Combine spectra
        res <- vector("list", len_z)
        hwsp <- hws + 1L
        for (i in seq_along(z)) {
            res[[i]] <- meanMzInts(z[windowIndices(i, hws, len_z)],
                                       weighted = wght, intensityFun = intF,
                                       main = hwsp - (i <= hws) * (hwsp - i),
                                       mzd = mzd, timeDomain = timeD, ppm = ppm,
                                       unionPeaks = FALSE)
        }
        res
    }, intF = intensityFun, wght = weighted, hws = as.integer(halfWindowSize),
    mzd = mzd, timeD = timeDomain, ppm = ppm, BPPARAM = BPPARAM)
    new_sp <- unsplit(new_sp, fromFile(x))
    names(new_sp) <- featureNames(x)
    x@assayData <- list2env(new_sp)
    lockEnvironment(x@assayData, bindings = TRUE)
    if (validObject(x))
        x
}

plotXIC_MSnExp <- function(x, ...) {
    ## Restrict to MS level 1
    x <- filterMsLevel(x, 1L)
    if (!length(x))
        stop("No MS1 data available")
    fns <- basename(fileNames(x))
    if (isMSnbaseVerbose())
        message("Retrieving data ...", appendLF = FALSE)
    f <- factor(fileNames(x), levels = fileNames(x))
    x <- splitByFile(x, f = force(f))
    x <- lapply(x, as, "data.frame")
    if (isMSnbaseVerbose())
        message("OK")
    ## Check if we are greedy and plot a too large area
    if (any(unlist(lapply(x, nrow)) > 20000))
        warning("The MS area to be plotted seems rather large. It is suggested",
                " to restrict the data first using 'filterRt' and 'filterMz'. ",
                "See also ?chromatogram and ?Chromatogram for more efficient ",
                "functions to plot a total ion chromatogram or base peak ",
                "chromatogram.",
                immediate = TRUE, call = FALSE)
    ## Define the layout.
    dots <- list(...)
    if (any(names(dots) == "layout")) {
        if (!is.null(dots$layout))
            layout(layout)
        dots$layout <- NULL
    } else
        layout(.vertical_sub_layout(length(x)))
    tmp <- mapply(x, fns, FUN = function(z, main, ...) {
        .plotXIC(x = z, main = main, layout = NULL, ...)
    }, MoreArgs = dots)
}
