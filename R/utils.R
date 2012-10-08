.get.amino.acids <- function() {
  get("amino.acids",envir=.MSnbaseEnv)
}

formatRt <- function(rt) {
  min <- floor(rt/60)
  sec <- round(rt-(min*60))
  return(paste(min,":",sec,sep=""))
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
    if (is.na(x[i]) | (x[i] == ""))
      x[i] <- x[i - 1]
  }
  return(x)
}

##' Return the name of variable \code{varname} in call \code{match_call}. 
##'
##' @title Return a variable name
##' @param match_call An object of class \code{call}, as returned by \code{match.call}. 
##' @param varname An \code{character} of length 1 which is looked up in \code{match_call}.
##' @return A \code{character} with the name of the variable passed as parameter
##' \code{varname} in parent close of \code{match_call}.
##' @examples
##' a <- 1
##' f <- function(x, y) 
##'  MSnbase:::getVariableName(match.call(), "x")
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
  if (n < 1)
    stop("'n' must be greater or equal than 1.")
  n <- min(n, nrow(X))
  X <- apply(X, 1, fun, ...)
  base::order(X, decreasing = TRUE)[1:n]
}

subsetBy <- function(X, groups, byIdx) {
  ans <- c()
  if ( is.null(dim(X)) || ncol(X) == 1 ) {
    X <- as.vector(X)
    for (l_i in unique(groups)) {
      X_i <- X[groups == l_i]
      j <- byIdx[[l_i]]
      ans <- c(ans, X_i[j])
    }    
  } else {
    for (l_i in unique(groups)) {
      X_i <- X[groups == l_i, ]
      j <- byIdx[[l_i]]
      ifelse(is.vector(X_i),
             ans <- base::rbind(ans, X_i),
             ans <- base::rbind(ans, X_i[j, ]))
    }
  }
  return(ans)
}



## Computes header from assay data by-passing cache
.header <- function(object) {
  if (length(object) == 0) 
    return(data.frame())
  if (all(msLevel(object) == 1)) {
    ln <- length(object)
    nas <- rep(NA, ln)
    hd <- list(file = fromFile(object),
               retention.time = rtime(object),
               precursor.mz = nas,
               precursor.intensity = nas,
               charge = nas,
               peaks.count = peaksCount(object),
               tic = tic(object),
               ionCount = ionCount(object),
               ms.level = msLevel(object),
               acquisition.number = acquisitionNum(object),
               collision.energy = nas)
  } else {
    ## tbl <- table(fromFile(object))
    ## idx <- as.numeric(unlist(apply(tbl, 1, function(x) 1:x)))
    hd <- list(file = fromFile(object),
               retention.time = rtime(object),
               precursor.mz = precursorMz(object),
               precursor.intensity = precursorIntensity(object),
               charge = precursorCharge(object),
               peaks.count = peaksCount(object),
               tic = tic(object),
               ionCount = ionCount(object),
               ms.level = msLevel(object),
               acquisition.number = acquisitionNum(object),
               collision.energy = collisionEnergy(object))
  }
  ## items are either a numeric or a list of integer() - keep former only
  sel <- sapply(hd, function(i) !is.list(i))
  hd <- as.data.frame(hd[sel]) 
  return(hd)
}

checkHeader <- function(object) {
  if (object@.cache$level == 0) {
    ans <- TRUE
  } else {
    cachedhd <- header(object)
    fromdata <- .header(object)
    ans <- identical(cachedhd, fromdata)
  }
  return(ans)
}


updateSpectrum2 <- function(x) {
  ## Version 0.2.0 of Spectrum has now a tic slot (MSnbase 1.5.3)
  newx <- new("Spectrum2")
  newx@merged <- x@merged
  newx@precScanNum <- x@precScanNum
  newx@precursorMz <-   x@precursorMz
  newx@precursorIntensity <- x@precursorIntensity
  newx@precursorCharge <- x@precursorCharge
  newx@collisionEnergy <- x@collisionEnergy   
  newx@msLevel <- x@msLevel
  newx@peaksCount <- x@peaksCount
  newx@rt <- x@rt
  newx@acquisitionNum <- x@acquisitionNum
  newx@scanIndex <- x@scanIndex
  newx@mz <- x@mz
  newx@intensity <- x@intensity
  newx@fromFile <- x@fromFile
  newx@centroided <- x@centroided
  newx@tic <- 0
  if (validObject(newx))
    return(newx)
}

updateMSnExp <- function(x) {
  for (fn in featureNames(x)) {
    .x <- x[[fn]]
    assign(fn, updateSpectrum2(.x), envir = assayData(x))
  }
  if (validObject(x))
    return(x)
}
 
