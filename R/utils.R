.get.amino.acids <- function() {
  get("amino.acids",envir=.MSnbaseEnv)
}

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
##' @seealso \code{\link{readMSData}}, \code{"\linkS4class{Spectrum1}"}, \code{"\linkS4class{Spectrum2}"}. 
##' @author Laurent Gatto
rawToSpectrum1 <- function(x,
                           i,
                           removePeaks=0,
                           clean=FALSE,
                           centroided=FALSE,
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
            scanIndex=x$scanindex,
            centroided=centroided)
  if (validObject(sp))
    return(sp)
}
  
rawToSpectrum2 <- function(x,
                           i,
                           removePeaks=0,
                           clean=FALSE,
                           centroided=FALSE,
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
            centroided=centroided,
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
  sapply(peakRanges,function(x) {
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
  sapply(zeroRanges,function(x){
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

makeImpuritiesMatrix <- function(x) {
  if (x==4) {
    M <- matrix(c(0.929,0.059,0.002,0.000,
                  0.020,0.923,0.056,0.001,
                  0.000,0.030,0.924,0.045,
                  0.000,0.001,0.040,0.923),
                nrow=4)
  } else {
    M <- diag(x)
  }
  colnames(M) <- paste("reporter",1:x,sep=".")
  rownames(M) <- paste("% reporter",1:x)
  corrfactors <- edit(M)
  invisible(corrfactors)
}

utils.removePrecMz <- function(spectrum, precMz=NULL,width=2) {
  ## Contributed by Guangchuang Yu for the plotMzDelta QC
  ## Additional modifications: setting peaks to 0 and clean argument
  if (is.null(precMz)) 
    precMz <- precursorMz(spectrum)
  if (!is.numeric(precMz)) 
    stop("precMz must either 'NULL' or numeric.")
  if (length(precMz) > 2) 
    stop ("precMz must a vector of length 1 or 2.")
  if (length(precMz) == 1) 
    precMz <- c(precMz-width, precMz+width)
  mz <- mz(spectrum)
  i <- intensity(spectrum)
  idx <- which(mz > precMz[1] & mz < precMz[2])
  spectrum@intensity[idx] <- 0
  return(spectrum)
}

utils.getMzDelta <- function(spectrum, percentage) {
  ## Computes the m/z differences between all the 
  ## 'percentage' top intensity peaks in a spectrum
  ## Contributed by Guangchuang Yu for the plotMzDelta QC
  mz <- mz(spectrum)
  i <- intensity(spectrum)
  idx <- order(i, decreasing=TRUE)
  tops <- idx[1:floor(length(idx) * percentage)] ## top 'percentage' of peaks
  mz.filtered <- mz[tops]
  delta <- c()   # mass delta
  while(length(mz.filtered) > 1) {
    m <- mz.filtered[1]
    mz.filtered <- mz.filtered[-1]
    d <- abs(mz.filtered-m)
    delta <- c(delta, d)
  }
  return(delta)
}	


fillUp <- function(x) {
  if (!any(is.na(x)) & !any(x != "")) 
    return(x)
  for (i in 2:length(x)) {
    if (is.na(x[i])) 
      x[i] <- x[i - 1]
    if (x[i] == "") 
      x[i] <- x[i - 1]
  }
  return(x)
}

##' Return the name of variable \code{varname} in call \code{match_call}. 
##'
##' @title Retirn variable name
##' @param match_call An object of class \code{call}, as returned by \code{match.call}. 
##' @param varname An \code{character} of length 1 which is looked up in \code{match_call}.
##' @return A \code{character} with the name of the variable passed as parameter
##' \code{varname} in parent close of \code{match_call}.
##' @examples
##' a <- 1
##' f <- function(x, y) 
##'  getVariableName(match.call(), "x")
##' f(x = a)
##' f(y = a)
##' @author Laurent Gatto
getVariableName <- function(match_call, varname) {
  match_call <- as.list(match_call)
  varname <- varname[1]
  mcx <- match_call[[varname]]
  while (any(sapply(mcx, length) != 1))
    mcx <- unlist(lapply(mcx, as.list))
  tail(as.character(mcx), n = 1)
}



##
## utils for topN method: getTopIdx and subsetBy
##

getTopIdx <- function(X, n, fun, ...) {
  ## Rows of X are first summerised using fun.
  ## Indices of the n highest values of vector X
  ## are then returned.
  ## input X: matrix [m,l]
  ## output: numeric of length min(n, nrow(x))  
  ## If (l == 1), fun does not have any effect.
  ## Otherwise, fun is required to keep the features
  ## grouped into rows.
  n <- min(n, nrow(X))
  X <- apply(X, 1, fun, ...)
  base::order(X, decreasing = TRUE)[1:n]
}

subsetBy <- function(X, groups, byIdx) {
  ans <- c()
  if ( !is.null(dim(X)) ) {
    for (l_i in unique(groups)) {
      X_i <- X[groups == l_i, ]
      j <- byIdx[[l_i]]
      ifelse(is.vector(X_i),
             ans <- base::rbind(ans, X_i),
             ans <- base::rbind(ans, X_i[j, ]))
    }
  } else {
    for (l_i in unique(groups)) {
      X_i <- X[groups == l_i]
      j <- byIdx[[l_i]]
      ans <- c(ans, X_i[j])
    }
  }
  return(ans)
}
