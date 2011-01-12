show.Spectrum <- function(spectrum) {
  if (msLevel(spectrum)==1) show.Spectrum1(spectrum)
  else show.Spectrum2(spectrum)
}
    

removePeaks.Spectrum <- function(spectrum,t="min") {
  if (t=="min") 
    t <- min(intensity(spectrum)[intensity(spectrum)>0])
  ints <- utils.removePeaks(spectrum@intensity,t)
  spectrum@intensity <- ints
  return(spectrum)
}


clean.Spectrum <- function(spectrum,updatePeaksCount=TRUE) {
  keep <- utils.clean(spectrum@intensity)
  spectrum@intensity <- spectrum@intensity[keep]
  spectrum@mz <- spectrum@mz[keep]
  spectrum@peaksCount <- length(spectrum@intensity)
  return(spectrum)
}

quantify.Spectrum <- function(spectrum,reporters,method) {
  ## Parameters:
  ##  spectrum: object of class Spectrum
  ##  reporters: object of class ReporterIons
  ## Return value:
  ##  a names list of length 2 with
  ##   peakArea: named numeric of length length(reporters)
  ##   curveStats: a length(reporters)x7 data frame 
  peakArea <- vector("numeric",length(reporters))
  names(peakArea) <- reporterNames(iTRAQ4)
  curveStats <- c()
  for (i in 1:length(reporters)) {
    ## Curve statistics
    dfr <- curveData(spectrum,reporters[i])
    ##  dfr:     mz int
    ##  1  114.1023   0
    ##  2  114.1063   2
    ##  3  114.1102   3
    ##  4  114.1142   4
    ##  ...
    maxInt <- max(dfr$int)
    nMaxInt <- sum(dfr$int==maxInt)
    baseLength <- nrow(dfr)
    mzRange <- range(dfr$mz)
    curveStats <- rbind(curveStats,
                        c(maxInt,nMaxInt,baseLength,
                          mzRange,
                          reporters[i]@mz,precursorMz(spectrum)))
    ## Quantification
    if (method=="trapezoidation") {
      ## Quantify reporter ions calculating the area
      ## under the curve by trapezoidation
      if (nrow(dfr)==1) {
        peakArea[i] <- 0
      } else {
        n <- nrow(dfr)
        x <- vector(mode = "numeric", length = n)
        for (j in 1:n) {
          k <- (j%%n) + 1
          x[j] <- dfr$mz[j] * dfr$int[k] - dfr$mz[k] * dfr$int[j]
        }
        peakArea[i] <- abs(sum(x)/2)
      }
    } 
    else if (method=="sum") {
      ## Quantify reporter ions using sum of data points of the peak
      peakArea[i] <- sum(dfr$int)
    }
    else if (method=="max") {
      ## Quantify reporter ions using max peak intensity
      peakArea[i] <- max(dfr$int)
    }
  }
  colnames(curveStats) <- c("maxInt","nMaxInt","baseLength",
                            "lowerMz","upperMz",
                            "reporter","precursor")
  curveStats <- as.data.frame(curveStats)  
  return(list(peakArea=peakArea,
              curveStats=curveStats))
}



## curveStats.Spectrum <- function(spectrum,reporters) {
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
    warning("Only returning data for first reporter ion")
    reporter <- reporter[1]
  }
  bp <- getCurveWidth(spectrum,reporter)
  if (any(is.na(bp))) {
    ## returning 'fake' data
    ## intensity 0 for 'missing' reporter
    return(data.frame(mz=reporter@mz,int=0))
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
  ## Warnings: the function returns warnings if the 
  ##  mz[indeces] range outside of the original window
  ##  reporter in the reporters object
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
      warning("[getCurveData] No data for for precursor ",spectrum@precursorMz," reporter ",m[i])
      xlwr[i] <- xupr[i] <- NA
    } else {
      ymax <- max(int[region])
      xmax <- which((int %in% ymax) & region)      
      xlwr[i] <- min(xmax) ## if several max peaks
      xupr[i] <- max(xmax) ## if several max peaks
      ylwr <- yupr <- ymax
      while (ylwr!=0) {
        xlwr[i] <- xlwr[i]-1
        ylwr <- int[xlwr[i]]
      }
      if (mz[xlwr[i]]<lwr[i])
        warning("Peak base for precursor ",spectrum@precursorMz,
                " reporter ",m[i],":\n   ",mz[xlwr[i]],"<",m[i],"-",
                reporters@width)
      while (yupr!=0) {
        xupr[i] <- xupr[i]+1
        yupr <- int[xupr[i]]
      }
      if (mz[xupr[i]]>upr[i])
        warning("Peak base for precursor ",spectrum@precursorMz,
                " reporter ",m[i],":\n   ",mz[xlwr[i]],">",m[i],"+",
                reporters@width)
      ## Updating xlwr and xupr [*]
      xlwr[i] <- xlwr[i]+1
      if (xupr[i]==length(int))
        xupr <- xupr-1
    }
  }
  return(list(lwr=xlwr,upr=xupr))
}

trimMz.Spectrum <- function(x,mzlim,updatePeaksCount=TRUE) {
  mzmin <- min(mzlim)
  mzmax <- max(mzlim)
  sel <- x@mz>mzmin & x@mz<mzmax
  if (sum(sel)==0) {
    warning(paste("No data points between ",mzmin," and ",mzmax,
                  " for precursor ",precursorMz(x),
                  ".\nLeaving data as is.",sep=""))
    return(x)
  }
  x@mz <- x@mz[sel]
  x@intensity <- x@intensity[sel]
  if (updatePeaksCount)
    x@peaksCount <- as.integer(length(x@intensity))
  return(x)
}


