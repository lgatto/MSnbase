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
            scanindex=x$scanindex[i],
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
  ## curve of max height 't' set t zero.
  ##
  ## Example:
  ## The following two curves will be removed
  ##   t - - - - + - - - - - -      
  ##           +  +  or  +++   
  ##   0 - - +    + - - +   + -
  ##
  ## Condition to remove peaks/curve:
  ## For a candidat having an intensity t at position j,
  ## this functions checks that intensities at previous
  ## (indices i below) and successive (indices k below)
  ## are monotonically decreasing 
  ## index notation:  i  <- j -> k
  ##
  ## Note:
  ## To avoid crashes in the two while loops in case
  ## int[1] or int[length(int)] are not 0, we set
  ## int <- c(0,int,0) and return int[2:(length(int)-1)].
  ## This assures that curves below 't' at the
  ## very beginning or end of int will also be removed.
  ##
  int <- c(0,int,0)
  candidats <- which(int<=t & int>0)
  for (j in candidats) {
    setToZero <- TRUE ## if curve meets condition
    i <- j - 1
    k <- j + 1
    ## test increasing monotonicity before candidate
    while(int[i]>0 & setToZero) {
      if (int[i]<=int[i+1]) {
        i <- i - 1
      } else {
        setToZero <- FALSE
      }
    }
    ## test decreasing monotonicity after candidate
    while(int[k]>0 & setToZero) {
      if (int[k]<=int[k-1]) {
        k <- k + 1
      } else {
        setToZero <- FALSE
      }
    }
    if (setToZero)
      int[i:k] <- 0
  }
  return(int[2:(length(int)-1)])
}  


utils.clean <- function(x) {
  ## Given a vector of numerics, this function
  ## returns a vector of logicals setting FALSE
  ## where zero is found in x when this zero is
  ## not adjacent to a non-zero value.
  ## input:  0 0 1 1 1 0 0 0 0 1 0 1 1 1 0 0
  ## output: F T T T T T F F T T T T T T T F
  ## See vignette for more details
  n <- length(x)
  b <- as.logical(rep(1,n)) ## initialise to TRUE
  if (x[1]==0 & x[2]==0)
    b[1] <- FALSE
  for (i in 2:(n-1)) {
    if (sum(x[(i-1):(i+1)])==0)
      b[i] <- FALSE
  }
  if (x[n-1]==0 & x[n]==0)
    b[n] <- FALSE
  return(b)
}

area <- function(pts,verbose=FALSE) {
  ## pl: 2x2 data frame or matrix
  ##  with points along the rows
  ##  and x and y coordinations in cols 1 and 2
  if (any(dim(pts)!=c(2,2)))
    stop("requires 2(points)x2(x,y coordinates) matrix or data frame as input")
  area <- 0
  x <- pts[,1]
  y <- pts[,2]
  ## triangle
  b <- abs(x[1]-x[2])
  h <- abs(y[1]-y[2])
  area <- area + ((b*h)/2)
  if (verbose)
    cat(" trianlge area:",area,"\n")
  if (min(y)!=0) { ## there is a rectangle
    h <- min(y)
    if (verbose)
      cat(" rectangle area:",b*h,"\n")
    area <- area + (h*b)
  }
  if (verbose)
    cat(" total area:",area,"\n")
  return(area)
}


zoom <- function(x,w=0.05) {
  new("ReporterIons",
      mz=x,
      width=w,
      col=rep("grep",length(x)))
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

pseudo3dplot <- function(hx) {
  require(ggplot2)
  p <- ggplot(hx,aes(retention.time,precursor.mz)) +
    geom_point(aes(color=peaks.count)) + 
      scale_colour_gradientn(colour=colorRampPalette(c("grey","blue","red","yellow"))(100), 
                             breaks=seq(min(hx$peaks.count),max(hx$peaks.count),length=8))
  return(p)
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
