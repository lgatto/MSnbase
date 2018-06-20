##' Objects of class `MSnExp` should be created with the `readMsData`
##' function. This simple constructor is for test data generation
##' purpose.
##'
##' @title Simple `MSnExp` constructor.
##' @param x A `list` of `Spectrum` objects.
##' @param fData A feature meta-data `data.frame`. Default is `NULL`.
##' @param pData A sample meta-data `data.frame`. Default is `NULL`.
##' @return An `MSnExp` object.
##' @author Laurent Gatto
MSnExp <- function(x, fData = NULL, pData = NULL) {
    stopifnot(inherits(x, "list"))
    x <- lapply(x, function(xx) {
        xx@fromFile <- 1L
        xx
    })
    if (is.null(names(x)))
        names(x) <- seq_len(length(x))
    if (is.null(fData))
        fData <- data.frame(row.names = names(x))
    if (is.null(pData))
        pData <- data.frame(row.names = "1")
    pd <- new("MSnProcess", files = NA_character_)
    e <- as.environment(x)
    new("MSnExp",
        assayData = e,
        phenoData = new("NAnnotatedDataFrame"),
        featureData = new("AnnotatedDataFrame", data = fData),
        processingData = pd)

}

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
#' consecutive scans are aggregated using the `intensityFun` and `mzFun`.
#' m/z values of intensities from consecutive scans will never be exactly
#' identical, even if they represent signal from the same ion. The function
#' determines thus internally a similarity threshold based on differences
#' between m/z values within and between spectra below which m/z values are
#' considered to derive from the same ion. For robustness reasons, this
#' threshold is estimated on the 100 spectra with the largest number of
#' m/z - intensity pairs (i.e. mass peaks).
#'
#' See [combineSpectra()] for details.
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
                                       mzFun = base::mean,
                                       intensityFun = base::mean,
                                       mzd = NULL,
                                       timeDomain = TRUE,
                                       BPPARAM = bpparam()){
    if (!is(x, "MSnExp"))
        stop("'x' has to be a 'MSnExp' or an 'OnDiskMSnExp'")
    if (is(x, "OnDiskMSnExp"))
        x <- as(x, "MSnExp")
    ## Combine spectra per file
    new_sp <- bplapply(split(spectra(x), fromFile(x)), FUN = function(z, intF,
                                                                      mzF, hws,
                                                                      mzd,
                                                                      timeD) {
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
            res[[i]] <- combineSpectra(z[windowIndices(i, hws, len_z)],
                                       mzFun = mzF, intensityFun = intF,
                                       main = hwsp - (i <= hws) * (hwsp - i),
                                       mzd = mzd, timeDomain = timeD)
        }
        res
    }, intF = intensityFun, mzF = mzFun, hws = as.integer(halfWindowSize),
    mzd = mzd, timeD = timeDomain, BPPARAM = BPPARAM)
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
    x <- as(x, "data.frame")
    x <- split(x, x$file)
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



}

##' @title Combine profile-mode spectra signals
##'
##' @description
##'
##' Combine signals from profile-mode spectra into the values for the
##' *main* spectrum. This can improve centroiding of profile-mode data
##' by increasing the signal-to-noise ratio.
##'
##' @details
##'
##' The m/z values of the same ion in consecutive scans (spectra) of a LCMS run
##' will not be identical. Assuming that this random variation is much smaller
##' than the resolution of the MS instrument (i.e. the difference between
##' m/z values within each single spectrum), m/z value groups are defined
##' across the spectra and those containing m/z values of the `main` spectrum
##' are retained. The maximum allowed difference between m/z values for the
##' same ion is estimated as in [estimateMzScattering()]. Alternatively it is
##' possible to define this maximal m/z difference with the `mzd` parameter.
##' All m/z values with a difference smaller than this value are combined to
##' a m/z group.
##' Intensities and m/z values falling within each of these m/z groups are
##' aggregated using the `intensity_fun` and `mz_fun`, respectively. It is
##' highly likely that all QTOF profile data is collected with a timing circuit
##' that collects data points with regular intervals of time that are then later
##' converted into m/z values based on the relationship `t = k * sqrt(m/z)`. The
##' m/z scale is thus non-linear and the m/z scattering (which is in fact caused
##' by small variations in the time circuit) will thus be different in the lower
##' and upper m/z scale. m/z-intensity pairs from consecutive scans to be
##' combined are therefore defined by default on the square root of the m/z
##' values. With `timeDomain = FALSE`, the actual m/z values will be used.
##'
##' @param object An instance of class `MSnExp`
##'
##' @param fcol `character(1)` indicating which feature variable to use
##'     to combine spectra. Will be coerced as a `factor`.
##'
##' @param mzFun `function` to aggregate the m/z values per m/z group. Should be
##'     a function or the name of a function. The function is expected to
##'     return a `numeric(1)`. For `mzFun = "weighted.mean"` (note
##'     that the *name* of the function is passed!) the resulting m/z is
##'     determined as an intensity-weighted mean of spectras' m/z values.
##'
##' @param intensityFun `function` to aggregate the intensity values per m/z
##'     group. Should be a function or the name of a function. The function is
##'     expected to return a `numeric(1)`.
##'
##' @param mzd `numeric(1)` defining the maximal m/z difference below which
##'     values are grouped. If not provided, this value is estimated from the
##'     distribution of differences of m/z values from the provided spectra
##'     (see details).
##'
##' @param timeDomain `logical(1)` whether definition of the m/z values to be
##'     combined into one m/z is performed on m/z values
##'     (`timeDomain = FALSE`) or on `sqrt(mz)` (`timeDomain = TRUE`).
##'     Profile data from TOF MS instruments should be aggregated based
##'     on the time domain (see details). Note that a pre-defined `mzd` should
##'     also be estimated on the square root of m/z values if
##'     `timeDomain = TRUE`.
##'
##' @return`combineSpectra` returns an new `MSnExp` object where the
##'     combined spectra m/z and intensity values representing the
##'     aggregated values across the provided spectra. The returned
##'     spectrum contains the same number of m/z and intensity pairs
##'     than the spectrum with index `main` in `x`, also all other
##'     related information is taken from this spectrum.
##'
##' @author Johannes Rainer, Sigurdur Smarason, Laurent Gatto
##'
##' @seealso
##'
##' [estimateMzScattering()] for a function to estimate m/z scattering
##' in consecutive scans.
##'
##' [estimateMzResolution()] for a function estimating the m/z resolution of
##' a spectrum.
##'
##' [combineSpectraMovingWindow()] for the function to combine consecutive
##' spectra of an `MSnExp` object using a moving window approach.
##'
##' @md
##'
##' @examples
##'
##' set.seed(123)
##' mzs <- seq(1, 20, 0.1)
##' ints1 <- abs(rnorm(length(mzs), 10))
##' ints1[11:20] <- c(15, 30, 90, 200, 500, 300, 100, 70, 40, 20) # add peak
##' ints2 <- abs(rnorm(length(mzs), 10))
##' ints2[11:20] <- c(15, 30, 60, 120, 300, 200, 90, 60, 30, 23)
##' ints3 <- abs(rnorm(length(mzs), 10))
##' ints3[11:20] <- c(13, 20, 50, 100, 200, 100, 80, 40, 30, 20)
##'
##' ## Create a testing MSnExp object with 3 spectra
##' sp1 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
##'     intensity = ints1, rt = 1, centroided = TRUE)
##' sp2 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.01),
##'     intensity = ints2, rt = 2, centroided = TRUE)
##' sp3 <- new("Spectrum1", mz = mzs + rnorm(length(mzs), sd = 0.009),
##'     intensity = ints3, rt = 3, centroided = TRUE)
##' x <- MSnExp(list(sp1 = sp1, sp2 = sp2, sp3 = sp3))
##' ## combine all spectra
##' fData(x)$all <- rep(1, 3)
##' x123 <- combineSpectra(x, "all")
##' ## combine two first spectra
##' fData(x)$sp12 <- c(1, 1, 2)
##' x12 <- combineSpectra(x, "sp12")
##' ## combine two last spectra
##' fData(x)$sp23 <- c(1, 2, 2)
##' x23 <- combineSpectra(x, "sp23")
##' ## combine first and last spectra
##' fData(x)$sp13 <- c(1, 2, 1)
##' x13 <- combineSpectra(x, "sp13")
##' ## None is combined
##' fData(x)$none <- c(1, 2, 3)
##' combineSpectra(x, "none")
##' ## Visualisation of results
combineSpectra <- function(object, fcol,
                           mzFun = base::mean,
                           intensityFun = base::mean,
                           timeDomain = TRUE) {
    stopifnot(inherits(object, "MSnExp"))
    stopifnot(fcol %in% fvarLabels(object))
    if (length(unique(msLevel(object))) != 1)
        stop("Can only combine spectra with the same MS level")
    fn0 <- featureNames(object)
    .by <- factor(fData(object)[, fcol], levels = unique(fData(object)[, fcol]))
    x <- spectra(object)
    res <- tapply(x, .by, FUN = .combineSpectra,
                  mzFun = mzFun,
                  intensityFun = intensityFun,
                  timeDomain = timeDomain)
    sel <- !duplicated(.by)
    object2 <- object[which(sel)]
    names(res) <- fn0[sel]
    new("MSnExp",
        assayData = as.environment(res),
        featureData = featureData(object2),
        phenoData = phenoData(object2),
        processingData = processingData(object2),
        experimentData = experimentData(object2),
        protocolData = protocolData(object))
}

##' @param x `list` of `Spectrum` objects.
##' @param main `integer(1)` defining the *main* spectrum, i.e. the spectrum
##'     which m/z and intensity values get replaced and is returned.
##' @return
##' @md
##' @rdname combineSpectra
##'
##' ## Combine the spectra
##' sp_agg <- .combineSpectra(list(sp1, sp2, sp3))
##'
##' ## Plot the spectra before and after combining
##' par(mfrow = c(2, 1), mar = c(4.3, 4, 1, 1))
##' plot(mz(sp1), intensity(sp1), xlim = range(mzs[5:25]), type = "h", col = "red")
##' points(mz(sp2), intensity(sp2), type = "h", col = "green")
##' points(mz(sp3), intensity(sp3), type = "h", col = "blue")
##' plot(mz(sp_agg), intensity(sp_agg), xlim = range(mzs[5:25]), type = "h", col = "black")
.combineSpectra <- function(x,
                            mzFun = base::mean,
                            intensityFun = base::mean,
                            main = floor(length(x) / 2L) + 1L,
                            mzd,
                            timeDomain = TRUE) {
    if (length(x) == 1)
        return(x[[1]])
    if (length(unique(unlist(lapply(x, function(z) z@msLevel)))) != 1)
        stop("Can only combine spectra with the same MS level")
    mzs <- lapply(x, function(z) z@mz)
    mzs_lens <- base::lengths(mzs)
    mzs <- unlist(mzs, use.names = FALSE)
    mz_order <- base::order(mzs)
    mzs <- mzs[mz_order]
    if (timeDomain)
        mz_groups <- .group_mz_values(sqrt(mzs), mzd = mzd)
    else
        mz_groups <- .group_mz_values(mzs, mzd = mzd)
    if (length(unique(mz_groups)) < length(x[[main]]@mz))
        stop("Got less m/z groups than m/z values in the original spectrum. ",
             "Most likely the data is not profile-mode LCMS data.")
    ## Want to keep only those groups with a m/z from the main spectrum.
    ## vectorized version from @sgibb
    is_in_main <- rep.int(seq_along(mzs_lens), mzs_lens)[mz_order] == main
    keep <- mz_groups %in% mz_groups[is_in_main]
    ## Keep only values for which a m/z in main is present.
    mz_groups <- mz_groups[keep]
    mzs <- mzs[keep]
    ints <- unlist(base::lapply(x, function(z) z@intensity))[mz_order][keep]
    ## Create result.
    new_sp <- x[[main]]
    ## Support also weighted.mean:
    if (is.character(mzFun) && mzFun == "weighted.mean") {
        intsp <- split(ints, mz_groups)
        new_sp@mz <- base::mapply(split(mzs, mz_groups), intsp,
                                  FUN = function(mz_vals, w)
                                      stats::weighted.mean(mz_vals, w + 1,
                                                           na.rm = TRUE),
                                  USE.NAMES = FALSE, SIMPLIFY = TRUE)
        new_sp@intensity <- base::vapply(intsp, FUN = intensityFun,
                                         FUN.VALUE = numeric(1),
                                         USE.NAMES = FALSE)
    } else {
        new_sp@mz <- base::vapply(split(mzs, mz_groups), FUN = mzFun,
                                  FUN.VALUE = numeric(1), USE.NAMES = FALSE)
        new_sp@intensity <- base::vapply(split(ints, mz_groups),
                                         FUN = intensityFun,
                                         FUN.VALUE = numeric(1),
                                         USE.NAMES = FALSE)
    }
    if (is.unsorted(new_sp@mz))
        stop("m/z values of combined spectrum are not ordered")
    new_sp@peaksCount <- length(new_sp@mz)
    new_sp
}
