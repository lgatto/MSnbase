removePeaks_Spectrum <- function(spectrum,t="min") {
  if (t=="min")
    t <- min(intensity(spectrum)[intensity(spectrum)>0])
  if (!is.numeric(t))
    stop("'t' must either be 'min' or numeric.")
  ints <- utils.removePeaks(spectrum@intensity,t)
  spectrum@intensity <- ints
  return(spectrum)
}


clean_Spectrum <- function(spectrum, all, updatePeaksCount = TRUE) {
  keep <- utils.clean(spectrum@intensity, all)
  spectrum@intensity <- spectrum@intensity[keep]
  spectrum@mz <- spectrum@mz[keep]
  spectrum@peaksCount <- length(spectrum@intensity)
  return(spectrum)
}

quantify_Spectrum <- function(spectrum, method,
                              reporters, strict) {
  ## Parameters:
  ##  spectrum: object of class Spectrum
  ##  method: methods for how to quantify the peaks
  ##  reporters: object of class ReporterIons
  ##  strict: use full width of peaks or a limited range around apex
  ## Return value:
  ##  a names list of length 2 with
  ##   peakQuant: named numeric of length length(reporters)
  ##   curveStats: a length(reporters) x 7 data frame
  peakQuant <- vector("numeric", length(reporters))
  names(peakQuant) <- reporterNames(reporters)
  curveStats <- c()
  for (i in 1:length(reporters)) {
      ## Curve statistics
      dfr <- curveData(spectrum, reporters[i])
      ##  dfr:     mz int
      ##  1  114.1023   0
      ##  2  114.1063   2
      ##  3  114.1102   3
      ##  4  114.1142   4
      ##  ...
      if (strict) {
          lowerMz <- round(mz(reporters[i]) - width(reporters[i]),3)
          upperMz <- round(mz(reporters[i]) + width(reporters[i]),3)
          selMz <- (dfr$mz >= lowerMz) & (dfr$mz <= upperMz)
          dfr <- dfr[selMz,]
          dfr <- rbind(c(min(dfr$mz),0),
                       dfr,
                       c(max(dfr$mz),0))
      }
      maxInt <- max(dfr$int)
      nMaxInt <- sum(dfr$int == maxInt)
      baseLength <- nrow(dfr)
      if (strict)
          baseLength <- baseLength-2
      mzRange <- range(dfr$mz)
      precMz <- precursorMz(spectrum)
      if (length(precMz) != 1)
          precMz <- NA
      curveStats <- rbind(curveStats,
                          c(maxInt,nMaxInt,baseLength,
                            mzRange,
                            reporters[i]@mz,precMz))
      ## Quantification
      if (method == "trapezoidation") {
          if (nrow(dfr) == 1) {
              if (!is.na(dfr$int))
                  warning(paste("Found only one mz value for precursor ",precursorMz(spectrum),
                                " and reporter ",reporterNames(reporters[i]),".\n",
                                "  If your data is centroided, quantify with 'max'.",sep=""))
          }
          ## Quantify reporter ions calculating the area
          ## under the curve by trapezoidation
          n <- nrow(dfr)
          ## - original code -
          ## x <- vector(mode = "numeric", length = n)
          ## for (j in 1:n) {
          ##     k <- (j%%n) + 1
          ##     x[j] <- dfr$mz[j] * dfr$int[k] - dfr$mz[k] * dfr$int[j]
          ## }
          ## peakQuant[i] <- abs(sum(x)/2)
          ## - updated -
          peakQuant[i] <- 0.5 * sum((dfr$mz[2:n] - dfr$mz[1:(n-1)]) *
                                    (dfr$int[2:n] + dfr$int[1:(n-1)]))
          ## area using zoo's rollmean, but seems slightly slower
          ## peakQuant[i] <- abs(sum(diff(dfr$int)*rollmean(dfr$mz,2)))
      } else if (method == "sum") {
          ## Quantify reporter ions using sum of data points of the peak
          peakQuant[i] <- sum(dfr$int)
      }  else if (method == "max") {
          ## Quantify reporter ions using max peak intensity
          peakQuant[i] <- max(dfr$int)
      }
  }
  colnames(curveStats) <- c("maxInt","nMaxInt","baseLength",
                            "lowerMz","upperMz",
                            "reporter","precursor")
  curveStats <- as.data.frame(curveStats)
  return(list(peakQuant=peakQuant,
              curveStats=curveStats))
}


## curveStats_Spectrum <- function(spectrum,reporters) {
##   curveStats <- c()
##   for (i in 1:length(reporters)) {
##     dfr <- curveData(spectrum,reporters[i])
##     maxInt <- max(dfr$int)
##     nMaxInt <- sum(dfr$int==maxInt)
##     baseLength <- nrow(dfr)
##     mzRange <- range(dfr$mz)
##     curveStats <- rbind(curveStats,
##                         c(maxInt,nMaxInt,baseLength,
##                           mzRange,
##                           reporters[i]@mz,precursorMz(spectrum)))
##   }
##   colnames(curveStats) <- c("maxInt","nMaxInt","baseLength",
##                             "lowerMz","upperMz",
##                             "reporter","precursor")
##   return(as.data.frame(curveStats))
## }

curveData <- function(spectrum,reporter) {
  ## Returns a data frame with mz and intensity
  ## values (as columns) for all the points (rows)
  ## in the reporter spectrum. The base of the
  ## curve is extracted by getCurveWidth
  ## Paramters:
  ##  spectrum: object of class Spectrum
  ##  reporter: object of class ReporterIons of length 1
  ## Return value:
  ##  a data frame
  ##           mz int
  ## 1  114.1023   0
  ## 2  114.1063   2
  ## 3  114.1102   3
  ## 4  114.1142   4
  ## ...
  if (length(reporter)!=1) {
    warning("[curveData] Only returning data for first reporter ion")
    reporter <- reporter[1]
  }
  bp <- getCurveWidth(spectrum,reporter)
  if (any(is.na(bp))) {
    return(data.frame(mz=reporter@mz,int=NA))
  } else {
    int <- intensity(spectrum)[bp$lwr[1]:bp$upr[1]]
    mz <- mz(spectrum)[bp$lwr[1]:bp$upr[1]]
    return(data.frame(cbind(mz,int)))
  }
}

getCurveWidth <- function(spectrum,reporters) {
  ## This function returns curve base indices
  ## from a spectrum object for all the reporter ions
  ## in the reporter object
  ## CHANGED IN VERSION 1.1.2
  ## NOT ANYMORE Warnings: the function returns warnings if the
  ## NOT ANYMORE  mz[indeces] range outside of the original window
  ## NOT ANYMORE  reporter in the reporters object
  ## Parameters:
  ##  spectrum: object of class Spectrum
  ##  reporters: object of class ReporterIons
  ## Return value:
  ##  list of length 2
  ##   - list$lwr of length(reporters) lower indices
  ##   - list$upr of length(reporters) upper indices
  m <- reporters@mz
  lwr <- m-reporters@width
  upr <- m+reporters@width
  mz <- spectrum@mz
  int <- spectrum@intensity
  ## if first/last int != 0, this function crashes in
  ## the while (ylwr!=0)/(yupr!=0) loops. Adding leading/ending data points
  ## to avoid this. Return values xlwr and xupr get updated accordingly [*].
  mz <- c(0,mz,0)
  int <- c(0,int,0)
  ## x... vectors of _indices_ of mz values
  ## y... intensity values
  xlwr <- xupr<- c()
  for (i in 1:length(m)) {
    region <- (mz>lwr[i] & mz<upr[i])
    if (sum(region,na.rm=TRUE)==0) {
      ## warning("[getCurveData] No data for for precursor ",spectrum@precursorMz," reporter ",m[i])
      xlwr[i] <- xupr[i] <- NA
    } else {
      ymax <- max(int[region])
      xmax <- which((int %in% ymax) & region)
      xlwr[i] <- min(xmax) ## if several max peaks
      xupr[i] <- max(xmax) ## if several max peaks
      if (!centroided(spectrum)) {
        ylwr <- yupr <- ymax
        while (ylwr!=0) {
          xlwr[i] <- xlwr[i]-1
          ylwr <- int[xlwr[i]]
        }
        ## if (mz[xlwr[i]]<lwr[i])
        ##   warning("Peak base for precursor ",spectrum@precursorMz,
        ##           " reporter ",m[i],": ",mz[xlwr[i]]," < ",m[i],"-",
        ##           reporters@width,sep="")
        while (yupr!=0) {
          xupr[i] <- xupr[i]+1
          yupr <- int[xupr[i]]
        }
        ## if (mz[xupr[i]]>upr[i])
        ##   warning("Peak base for precursor ",spectrum@precursorMz,
        ##           " reporter ",m[i],": ",mz[xlwr[i]]," > ",m[i],"+",
        ##           reporters@width,sep="")
      }
      ##
      ## --|--|--|--|--|--|--
      ##      1  2  3  4      indeces before c(0,mz,0)
      ##   1  2  3  4  5  6   indeces after  c(0,mz,0)
      ## lower index: if x_after > 1 x <- x -1 (else x <- x_after)
      ## upper index: always decrement by 1
      ##              if we have reached the last index (the 0), decrement by 2
      ##
      ## Updating xlwr, unless we reached the artificial leading 0
      if (xlwr[i]>1)
        xlwr[i] <- xlwr[i]-1
      ## Always updating xupr [*]
      if (xupr[i]==length(mz))
        xupr[i] <- xupr[i]-2
      xupr[i] <- xupr[i]-1
    }
  }
  return(list(lwr=xlwr,upr=xupr))
}

trimMz_Spectrum <- function(x,mzlim,updatePeaksCount=TRUE) {
  mzmin <- min(mzlim)
  mzmax <- max(mzlim)
  sel <- (x@mz >= mzmin) & (x@mz <= mzmax)
  if (sum(sel)==0) {
    warning(paste("No data points between ", mzmin, " and ", mzmax,
                  " for spectrum with acquisition number ", acquisitionNum(x),
                  ".\n Leaving data as is.", sep=""))
    return(x)
  }
  x@mz <- x@mz[sel]
  x@intensity <- x@intensity[sel]
  if (updatePeaksCount)
    x@peaksCount <- as.integer(length(x@intensity))
  return(x)
}

normalise_Spectrum <- function(object, method, value) {
  ints <- intensity(object)
  switch(method,
         max = div <- max(ints),
         sum = div <- sum(ints),
         value = div <- value)
  normInts <- ints/div
  object@intensity <- normInts
  if (validObject(object))
    return(object)
}

bin_Spectrum <- function(object, binSize=1L,
                         breaks=seq(floor(min(mz(object))),
                                    ceiling(max(mz(object))), by=binSize),
                         fun=sum) {
  fun <- match.fun(fun)
  nb <- length(breaks)
  n <- peaksCount(object)

  idx <- findInterval(mz(object), breaks)

  idx[which(idx < 1L)] <- 1L
  idx[which(idx > n)] <- n

  intensity <- double(length(breaks))
  intensity[unique(idx)] <- unlist(lapply(split(intensity(object), idx), fun))

  mz <- c((breaks[-nb]+breaks[-1L])/2L, breaks[nb])

  object@mz <- mz
  object@intensity <- intensity
  object@tic <- sum(intensity)
  object@peaksCount <- nb
  return(object)
}

bin_Spectra <- function(object1, object2, binSize=1L,
                        breaks=seq(floor(min(c(mz(object1), mz(object2)))),
                                   ceiling(max(c(mz(object1), mz(object2)))),
                                   by=binSize)) {
  return(list(bin_Spectrum(object1, binSize=binSize, breaks=breaks),
              bin_Spectrum(object2, binSize=binSize, breaks=breaks)))
}

#' calculate similarity between spectra (between their intensity profile)
#' @param x spectrum1 (MSnbase::Spectrum)
#' @param y spectrum2 (MSnbase::Spectrum)
#' @param fun similarity function (must take two spectra and ... as arguments)
#' @param ... further arguments passed to "fun"
#' @return double, similarity score
compare_Spectra <- function(x, y, fun=c("common", "cor", "dotproduct"), ...) {
  if (is.character(fun)) {
    fun <- match.arg(fun)
    if (fun == "cor" || fun == "dotproduct") {
      binnedSpectra <- bin_Spectra(x, y, ...)
      inten <- lapply(binnedSpectra, intensity)
      return(do.call(fun, inten))
    } else if (fun == "common") {
      return(numberOfCommonPeaks(x, y, ...))
    }
  } else if (is.function(fun)) {
    return(fun(x, y, ...))
  }
  return(NA)
}

#' plot spectrum1 vs spectrum2
#' @param spectra list, 2 MSnbase::Spectrum2 objects
#' @param sequences list, 2 character vectors containing the peptide sequence
#' [not used yet]
#' @param fragments list, 2 character vectors containing the fragment strings
#' [not used yet]
#' @param norm normalize?
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param xlim limits for x-axis
#' @param ylim limits for y-axis
#' @param tolerance double, allowed deviation to be considered as common peak
#' @param relative relative (or absolute) deviation
#' @param fragments.cex cex for fragments
#' @param legend.cex cex for legend
plot_Spectra <- function(spectra,
                         sequences,
                         fragments,
                         norm=TRUE,
                         xlab="mz", ylab="intensity",
                         xlim, ylim,
                         tolerance=25e-6,
                         relative=TRUE,
                         fragments.cex=0.5,
                         legend.cex=1, ...) {
  centroided <- sapply(spectra, centroided)

  if (all(!centroided)) {
    message("Your spectra are not centroided.")
  } else if (any(!centroided)) {
    warning("Only one spectrum is not centroided!")
  }

  if (norm) {
    spectra <- lapply(spectra, normalize)
  }
  if (missing(xlim)) {
    mass <- unlist(lapply(spectra, mz))
    xlim <- c(min(c(mass, 0)), max(c(mass, 0)))
  }

  if (missing(ylim)) {
    inten <- unlist(lapply(spectra, intensity))
    maxInten <- max(c(inten, 0))
    ylim <- c(-maxInten, maxInten)
  }

  common <- lapply(list(c(1, 2), c(2, 1)), function(x) {
    commonPeaks(spectra[[x[1]]], spectra[[x[2]]],
                method="highest", tolerance=tolerance, relative=relative)
  })

  plot(NA, type="h", col=1,
       xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim,
       ...)
  abline(h=0, col="#808080")

  orientation <- c(1, -1)
  text.pos <- c(3, 1)
  legend.pos <- c("topleft", "bottomleft")
  legend.prefix <- c("ident", "quant")
  ## colors: ColorBrewer RdYlBu c(9, 11, 3, 1)
  cols <- c("#74ADD1", "#313695", "#F46D43", "#A50026")
  pch <- c(NA, 19)

  for (i in seq(along=spectra)) {
    lines(mz(spectra[[i]]), orientation[i]*intensity(spectra[[i]]),
          type="h", col=cols[(i-1)*2+common[[i]]+1L], lwd=1.5)
    points(mz(spectra[[i]]), orientation[i]*intensity(spectra[[i]]),
           col=cols[(i-1)*2+common[[i]]+1L], pch=pch[common[[i]]+1L],
           cex=0.5)

    if (!missing(fragments) && length(fragments[[i]])) {
      text(mz(spectra[[i]]), orientation[i]*intensity(spectra[[i]]),
           fragments[[i]], pos=text.pos[i], offset=0.25,
           cex=fragments.cex, col="#808080")
    }

    label <- paste0("prec scan: ", precScanNum(spectra[[i]]))

    if (peaksCount(spectra[[i]])) {
      label <- paste0(label, ", prec mass: ", round(precursorMz(spectra[[i]]), 3),
                             ", prec z: ", precursorCharge(spectra[[i]]),
                             ", # common: ", sum(common[[i]]))
      if (!missing(sequences)) {
        label <- paste0(label, ", seq: ", sequences[[i]])
      }
    }

    legend(legend.pos[i], legend=label, bty="n", cex=legend.cex)

pickPeaks_Spectrum <- function(object, halfWindowSize = 2L,
                               method = c("MAD", "SuperSmoother"),
                               SNR = 0L, ...) {

  if (!peaksCount(object)) {
    warning("Your spectrum is empty. Nothing to pick.")
    return(object)
  }

  if (length(object@centroided) && object@centroided) {
    warning("Your spectrum is already centroided.")
    return(object)
  }

  ## estimate noise
  noise <- MALDIquant:::.estimateNoise(mz(object), intensity(object),
                                       method = match.arg(method), ...)

  ## find local maxima
  isLocalMaxima <- MALDIquant:::.localMaxima(intensity(object),
                                             halfWindowSize = halfWindowSize)

  ## include only local maxima which are above the noise
  isAboveNoise <- object@intensity > (SNR * noise[, 2L])

  peakIdx <- which(isAboveNoise & isLocalMaxima)

  object@mz <- object@mz[peakIdx]
  object@intensity <- object@intensity[peakIdx]
  object@peaksCount <- length(peakIdx)
  object@centroided <- TRUE

  if (validObject(object)) {
    return(object)
  }
}

smooth_Spectrum <- function(object,
                            method = c("SavitzkyGolay", "MovingAverage"),
                            halfWindowSize = 2L, ...) {

  if (!peaksCount(object)) {
    warning("Your spectrum is empty. Nothing to change.")
    return(object)
  }

  method <- match.arg(method)

  switch(method,
         "SavitzkyGolay" = {
           fun <- MALDIquant:::.savitzkyGolay
         },
         "MovingAverage" = {
           fun <- MALDIquant:::.movingAverage
         })

  object@intensity <- fun(object@intensity, halfWindowSize = halfWindowSize, ...)

  isBelowZero <- object@intensity < 0

  if (any(isBelowZero)) {
    warning("Negative intensities generated. Replaced by zeros.")
    object@intensity[isBelowZero] <- 0
  }

  if (validObject(object)) {
    return(object)
  }
}

