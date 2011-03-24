##' Converts seconds to min:sec format
##'
##' This function is used to convert retention times,
##' expressed in seconds, in the more human friendly
##' format mm:sec.
##'
##' @title Format Retention Time
##' @usage formatRt(rt)
##' @param rt retention in in seconds
##' @return a character string in mm:ss format 
##' @examples
##' formatRt(1524)
##' @author Laurent Gatto
##' @export 
##' @keywords chron utilities
formatRt <- function(rt) {
  min <- floor(rt/60)
  sec <- round(rt-(min*60))
  return(paste(min,":",sec,sep=""))
}

##' Convert spectrum data to an Spectrum object
##'
##' This function converts the i'th spectrum of a list produced
##' by \code{xcms:::rampRawData} (for \code{rawToSpectrum1}) or
##' \code{xcms:::rampRawDataMSn} (for \code{rawToSpectrum2}) to
##' object of class \code{"\linkS4class{Spectrum1}"} or
##' \code{"\linkS4class{Spectrum1}"}, respectively.
##' Rudimentary data processing like low intensity peaks and
##' zero intensity data points removal can optionally already be
##' applied at this stage.
##' 
##' @title Raw spectrum data to Spectrum object 
##' @aliases rawToSpectrum2
##' @usage
##' rawToSpectrum1(x,i,removePeaks,clean,updatePeaksCount,verbose)
##' rawToSpectrum2(x,i,removePeaks,clean,updatePeaksCount,verbose)
##' @param x a list as returned by \code{xcms:::rampRawData} (for \code{rawToSpectrum1})
##' or \code{xcms:::rampRawDataMSn} (for \code{rawToSpectrum2}).
##' @param i a numeric value (integer) indicating which spectrum to extract.
##' @param removePeaks a numeric value to define low-intensity peaks that should
##' be removed; indicate 0 (default) to keep data as is.
##' @param clean a logical value specifying whether MZ data points with 0 intensity
##' (if for instance removePeaks > 1) should be completely removed.
##' @param updatePeaksCount a logical value indicating if peak counts should
##' be updated (default is TRUE); should be set to TRUE.
##' @param verbose a logical value indicating if ouput is written when removing
##' low intensity peaks and cleaning data (default is FALSE).
##' @return an object of class \code{"\linkS4class{Spectrum1}"} for \code{rawToSpectrum1}
##' or \code{"\linkS4class{Spectrum1}"} for \code{rawToSpectrum2}.
##' @seealso \code{\link{readMzXMLData}}, \code{"\linkS4class{Spectrum1}"}, \code{"\linkS4class{Spectrum2}"}. 
##' @author Laurent Gatto
rawToSpectrum1 <- function(x,
                           i,
                           removePeaks=0,
                           clean=FALSE,
                           updatePeaksCount=TRUE,
                           verbose=verbose) {
  from <-x$scanindex[i]+1
  to <- x$scanindex[i]+x$tic[i]
  int <- x$intensity[from:to]
  mz <- x$mz[from:to]
  pc <- x$tic[i]
  if (removePeaks>0) {
    if (verbose) cat(" [rawToSpectru1] removing peaks <=",removePeaks,"\n")
    int <- utils.removePeaks(int,removePeaks)
    clean <- TRUE
  }
  if (clean) {
    if (verbose) cat(" [rawToSpectrum1] cleaning...\n")
    keep <- utils.clean(int)
    int <- int[keep]
    mz <- mz[keep]
    if (updatePeaksCount)
      pc <- length(int)
  }
  sp <- new("Spectrum1",
            mz=mz,
            intensity=int,
            peaksCount=as.integer(pc),
            polarity=x$polarity[i],
            rt=x$rt[i],
            acquisitionNum=x$acquisitionNum[i],
            scanIndex=x$scanindex)
  if (validObject(sp))
    return(sp)
}
  
rawToSpectrum2 <- function(x,
                           i,
                           removePeaks=0,
                           clean=FALSE,
                           updatePeaksCount=TRUE,
                           verbose=verbose) {
  from <-x$scanindex[i]+1
  to <- x$scanindex[i]+x$peaksCount[i]
  int <- x$intensity[from:to]
  mz <- x$mz[from:to]
  pc <- x$peaksCount[i]
  precursorMZ <- x$precursorMZ[i]
  if (removePeaks>0) {
    if (verbose) cat(" [rawToSpectru2] removing peaks <=",removePeaks,"\n")
    int <- utils.removePeaks(int,removePeaks)
    clean <- TRUE
  }
  pc <- length(int) ## peaksCount before cleaning
  if (clean) {
    if (verbose) cat(" [rawToSpectrum2] cleaning...\n")
    keep <- utils.clean(int)
    int <- int[keep]
    mz <- mz[keep]
    if (updatePeaksCount)
      pc <- length(int)
  }
  sp <- new("Spectrum2",
            rt=x$rt[i],
            acquisitionNum=x$acquisitionNum[i],
            precursorMz=x$precursorMZ[i],
            precursorIntensity=x$precursorIntensity[i],
            precursorCharge=x$precursorCharge[i],
            scanIndex=x$scanindex[i],
            collisionEnergy=x$collisionEnergy[i],
            peaksCount=as.integer(pc),
            mz=mz,
            intensity=int)
  if (validObject(sp))
    return(sp)
}

utils.removePeaks <- function(int,t) {
  ## Description:
  ## Given a vector of intensities 'int' and a threshold 't',
  ## this function returns vector of same length with all
  ## peaks of max height 't' set t zero.
  ## Example:
  ## The following three curves will be removed
  ##   t - - - - + - - - - - - - - + - + -     
  ##           +  +  or  +++     ++ +++ +
  ##   0 - - +    + - - +   + - + - - - +  
  ##  
  peakRanges <- IRanges(sapply(int,">",0))
  IRanges::sapply(peakRanges,function(x) {
    ## we get the indices of every peak in int
    if(all(int[x]<=t)) 
      int[x] <<- 0
  })
  return(int)
}

utils.clean <- function(x) {
  ## Given an numeric x, this function
  ## returns a logical b of length(x) where
  ## non-zero values and their direct
  ## 0s are TRUE so that x[b] has only
  ## non-zero values surrounded by it's
  ## original direct zero neighbours.
  ## Example
  ## x: 1 0 0 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0
  ## b: T T F T T T T T T T T T F T T T F F
  ##
  ## x[b]:     1 0 1 1 1 0 1 1 0 1 0
  n <- length(x)
  b <- as.logical(rep(1,n)) ## initialise to TRUE
  zeroRanges <- IRanges(sapply(x,"==",0))
  IRanges::sapply(zeroRanges,function(x){
    if (length(x)>2)
      b[x[2:(length(x)-1)]] <<- FALSE 
  })
  return(b)
}

zoom <- function(x,w=0.05) {
  new("ReporterIons",
      mz=x,
      width=w,
      name="xlim",
      reporterNames=paste("xlim",x,sep="."),
      col=rep("grey",length(x)))
}


getRatios <- function(x,log=FALSE) {
  ## x: a vector of numerics
  ## returns a vector of all xi/xj ratios
  x <- as.numeric(x)
  cmb <- combn(length(x),2)
  r <- numeric(ncol(cmb))
  for (i in 1:ncol(cmb)) {
    j <- cmb[1,i]
    k <- cmb[2,i]
    ifelse(log,r[i] <- x[j]-x[k],r[i] <- x[j]/x[k])
  }
  return(r)
}


getBins <- function(x) {
  bins <- numeric(length(x))
  bins[1] <- 1
  for (i in 2:length(x)) {
    ifelse(x[i]==x[i-1]+1,
           bins[i] <- bins[i-1],
           bins[i] <- bins[i-1]+1)
  }
  return(bins)
}

utils.plot2d <- function(object,z=c("tic","file","peaks.count","charge")) {
  z <- match.arg(z)
  stopifnot(c("retention.time","precursor.mz",z) %in% names(object))
  peaks.count <- charge <- retention.time <- precursor.mz <- NULL # to satisfy codetools
  p <- ggplot(object,aes(retention.time,precursor.mz)) + labs(colour=z)
  switch(z,
         tic = p <- p+ geom_point(aes(colour=tic)),
         peaks.count = p <- p+ geom_point(aes(colour=peaks.count)),
         file = p <- p+ geom_point(aes(colour=as.factor(file))),
         charge = p <- p+ geom_point(aes(colour=as.factor(charge))))
  print(p)
  invisible(p)
}
