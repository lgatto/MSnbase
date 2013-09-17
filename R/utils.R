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

utils.clean <- function(x, all) {
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
  if (all) {    
    b[x == 0] <- FALSE
  } else {
    zeroRanges <- IRanges(sapply(x,"==",0))

    sapply(zeroRanges,function(x){
      if (length(x)>2)
        b[x[2:(length(x)-1)]] <<- FALSE 
    })
  }
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

makeImpuritiesMatrix <- function(x, filename, edit = TRUE) {
  if (!missing(filename)) {
    m <- read.csv(filename, row.names = 1)
    x <- ncol(m)
    if (ncol(m) != nrow(m))
      stop(paste0("Problem reading impurity matrix. Not square.\n",
                  "Please read '?purityCorrect' for details."))
    ncharge <- x/2
    a <- (x/2)
    b <- (x/2) + 1
    res <- matrix(0, x, x)
    diag(res) <- 100 - rowSums(m)
    for (k in 1:ncharge) {
      diag(res[(1+k):x, 1:(x-k)]) <- m[(1+k):x, (a-k+1)]
      diag(res[1:(x-k), (1+k):x]) <- m[1:(x-k), (b+k-1)]
    }
    ## test <- matrix(0, 6, 6)
    ## diag(test) <- 100 - rowSums(m)
    ## diag(test[4:6, 1:3]) <- m[4:6, 1] ## col1: -3
    ## diag(test[3:6, 1:4]) <- m[3:6, 2] ## col2: -2 
    ## diag(test[2:6, 1:5]) <- m[2:6, 3] ## col3: -1
    ## diag(test[1:5, 2:6]) <- m[1:5, 4] ## col4: +1
    ## diag(test[1:4, 3:6]) <- m[1:4, 5] ## col5: +2
    ## diag(test[1:3, 4:6]) <- m[1:3, 6] ## col6: +3
    ## test <- test/100
    M <- res/100
  } else {
    if (x==4) {
      M <- matrix(c(0.929,0.059,0.002,0.000,
                    0.020,0.923,0.056,0.001,
                    0.000,0.030,0.924,0.045,
                    0.000,0.001,0.040,0.923),
                  nrow=4, byrow = TRUE)
    } else if (x == 6) {
      M <- matrix(c(0.939, 0.061, 0.000, 0.000, 0.000, 0.000,
                    0.005, 0.928, 0.067, 0.000, 0.000, 0.000,
                    0.000, 0.011, 0.947, 0.042, 0.000, 0.000,
                    0.000, 0.000, 0.017, 0.942, 0.041, 0.000,
                    0.000, 0.000, 0.000, 0.016, 0.963, 0.021,
                    0.000, 0.000, 0.000, 0.002, 0.032, 0.938),
                  nrow = 6, byrow = TRUE)
    } else {
      M <- diag(x)
    }
  }  
  colnames(M) <- paste("reporter", 1:x, sep=".")
  rownames(M) <- paste("% reporter", 1:x)
  if (edit)
    M <- edit(M)
  return(M)
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
 

cramer4 <- function(object, imp) {
  ## see Shadford et al. 2005, BMC Genomics
  if (missing(imp)) {
    impM <- matrix(c(0.0, 1.0, 5.9, 0.2,
                     0.0, 2.0, 5.6, 0.1,
                     0.0, 3.0, 4.5, 0.1,
                     0.1, 4.0, 3.5, 0.1),
                   nrow = 4, byrow = TRUE)
    colnames(impM) <- c("-2", "-1", "+1", "+2")
    rownames(impM) <- 114:117
    imp <- as.numeric(impM)
    names(imp) <- letters[1:length(imp)]
  }

  w <- (100 - (imp["a"] + imp["e"] + imp["i"] + imp["m"]))
  x <- (100 - (imp["b"] + imp["f"] + imp["j"] + imp["n"]))
  y <- (100 - (imp["c"] + imp["g"] + imp["k"] + imp["o"]))
  z <- (100 - (imp["d"] + imp["h"] + imp["l"] + imp["p"]))

  C <- matrix(c(w, imp["f"], imp["c"], 0,
                imp["i"], x, imp["g"], imp["d"],
                imp["m"], imp["j"], y, imp["h"],
                0, imp["n"], imp["k"], z),
              ncol = 4, byrow = TRUE)

  if (det(C) == 0) {
    warning("Determinant of C is 0, correction impossible")
    object@processingData@processing <-
      c(object@processingData@processing,
        "No impurity correction possible, det(C) is 0")    
  } else {    
    e <- exprs(object)
    res <- apply(e, 1, function(.e) {
      d1 <- matrix(c(.e,
                     imp["f"], x, imp["j"], imp["n"],
                     imp["c"], imp["g"], y, imp["k"],
                     0, imp["d"], imp["h"], z),
                   ncol = 4, byrow = TRUE)
      d2 <- matrix(c(w, imp["i"], imp["m"], 0,
                     .e,
                     imp["c"], imp["g"], y, imp["k"],
                     0, imp["d"], imp["h"], z),
                   ncol = 4, byrow = TRUE)
      d3 <- matrix(c(w, imp["i"], imp["m"], 0,
                     imp["f"], x, imp["j"], imp["n"],               
                     .e,
                     0, imp["d"], imp["h"], z),
                   ncol = 4, byrow = TRUE)
      d4 <- matrix(c(w, imp["i"], imp["m"], 0,
                     imp["f"], x, imp["j"], imp["n"],               
                     imp["c"], imp["j"], y, imp["k"],
                     .e),
                   ncol = 4, byrow = TRUE)
      res <- c(det(d1)/det(C),
               det(d2)/det(C),
               det(d3)/det(C),
               det(d4)/det(C))
      return(res)
    })
    res <- t(res)
    rownames(res) <- featureNames(object)
    colnames(res) <- sampleNames(object)
    object@processingData@processing <-
      c(object@processingData@processing,
        "Impurity correction using Cramer's rule.")
    exprs(object) <- res
  }
  if (validObject(object))
    return(object)
}

cramer6 <- function(x, imp) {
  if (missing(imp)) {
    imp <- c(0, 0, 0, 6.1, 0, 0,
             0, 0, 0.5, 6.7, 0, 0, 
             0, 0, 1.1, 4.2, 0, 0, 
             0, 0, 1.7, 4.1, 0, 0, 
             0, 0, 1.6, 2.1, 0, 0, 
             0, 0.2, 3.2, 2.8, 0, 0)
    names(imp) <- letters[1:length(imp)]
    impM <- matrix(imp, nrow = 6, byrow = TRUE)
    colnames(impM) <- c("-3", "-2", "-1", "+1", "+2", "+3")
    rownames(impM) <- 126:131
    imp <- as.numeric(imp)
  }
  return(FALSE)
}

     
getColsFromPattern <- function(x, pattern) {
  if (missing(pattern))
    stop("Pattern must not be missing.")
  if (nchar(pattern) != ncol(x))
    stop("The pattern must be equal to the number of columns.")
  pattern <- strsplit(pattern, "")[[1]]
  if (!all(unique(pattern) %in% c("0", "1")))
    stop("Pattern must be composed of '0' or '1' defining columns with or without 'NA's.")
  return(pattern == "1")  
}

getRowsFromPattern <- function(x, pattern) {
  cols <- getColsFromPattern(x, pattern)
  x2 <- x[, cols]
  apply(x2, 1, function(xx) !any(is.na(xx)))
}


nologging <- function(object, n = 1) {
  ## removes the last n entries from 
  ## object@processingData@processing
  l <- length(object@processingData@processing)
  x <- seq(l, length = n, by = -1)
  object@processingData@processing <-
    object@processingData@processing[-x]
  stopifnot(length(object@processingData@processing) == (l - n))
  if (validObject(object))
    return(object)
}

logging <- function(object, msg, date. = TRUE) {
  if (date.)
    msg <- paste0(msg, ": ", date())
  object@processingData@processing <-
    c(object@processingData@processing, msg)
  if (validObject(object))
    return(object)
}

##' Given a text spread sheet \code{f} and a \code{pattern} to
##' be matched to its header (first line in the file), the function
##' returns the matching columns names or indices of the
##' corresponding \code{data.frame}.
##' 
##' The function starts by reading the first line of the file (or connection)
##' \code{f} with \code{\link{readLines}}, then splits it
##' according to the optional \code{...} arguments (it is important to
##' correctly specify \code{\link{strsplit}}'s \code{split} character vector here)
##' and then matches \code{pattern} to the individual column names using
##' \code{\link{grep}}.
##'
##' Similarly, \code{getEcols} can be used to explore the column names and
##' decide for the appropriate \code{pattern} value.
##' 
##' These functions are useful to check the parameters to be provided to
##' \code{\link{readMSnSet2}}.
##' 
##' @title Returns the matching column names of indices.
##' @param f A connection object or a \code{character} string to
##' be read in with \code{readLines(f, n = 1)}. 
##' @param pattern A \code{character} string containing a regular expression
##' to be matched to the file's header.
##' @param ... Additional parameters passed to \code{\link{strsplit}}
##' to split the file header into individual column names.
##' @return Depending on \code{value}, the matching column names of
##' indices. In case of \code{getEcols}, a \code{character} of
##' column names.
##' @seealso \code{\link{readMSnSet2}}
##' @author Laurent Gatto
grepEcols <- function(f, pattern, ...) 
    grep(pattern, strsplit(readLines(f, 1), ...)[[1]])


##' @rdname grepEcols
getEcols <- function(f, ...)
    strsplit(readLines(f, 1), ...)[[1]]

MSnExp.size <- function(x) 
    object.size(x) + sum(unlist(unname(eapply(assayData(x),
                                              object.size)))) 

