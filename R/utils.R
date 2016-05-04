##' Returns a \code{data.frame} of amino acid properties: \code{AA},
##' \code{ResidueMass}, \code{Abbrev3}, \code{ImmoniumIonMass},
##' \code{Name}, \code{Hydrophobicity}, \code{Hydrophilicity},
##' \code{SideChainMass}, \code{pK1}, \code{pK2} and \code{pI}.
##'
##' @title Amino acids
##' @return A \code{data.frame}
##' @author Laurent Gatto
##' @examples
##' get.amino.acids()
get.amino.acids <- function()
    .get.amino.acids()

.get.amino.acids <- function() {
  get("amino.acids",envir=.MSnbaseEnv)
}

##' Returns a \code{double} of used atomic mass.
##'
##' @title Atomic mass.
##' @return A named \code{double}.
##' @author Sebastian Gibb
##' @examples
##' get.atomic.mass()
get.atomic.mass <- function()
    .get.atomic.mass()

.get.atomic.mass <- function() {
  get("atomic.mass",envir=.MSnbaseEnv)
}

formatRt <- function(rt) {
    ans <- NA
    if (is.numeric(rt)) {
        min <- floor(rt/60)
        sec <- round(rt-(min*60))
        ans <- paste(min,":",sec,sep="")
    } else if (is.character(rt)) {
        ans <- strsplit(rt, ":")
        ans <- sapply(ans, function(x) {
            x <- as.numeric(x)
            60 * x[1] + x[2]
        })
    } else {
        warning("Input must be numeric of character.")
    }
    return(ans)
}

utils.removePeaks_centroided <- function(int, t) {
    rmi <- int <= t
    int[rmi] <- 0
    int
}

utils.removePeaks <- function(int, t) {
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



utils.removePrecMz <- function(spectrum, precMz = NULL, width = 2) {
  ## Contributed by Guangchuang Yu for the plotMzDelta QC
  ## Additional modifications: setting peaks to 0 and clean argument
  if (is.null(precMz))
    precMz <- precursorMz(spectrum)
  if (!is.numeric(precMz))
    stop("precMz must either 'NULL' or numeric.")
  if (length(precMz) > 2)
    stop ("precMz must a vector of length 1 or 2.")
  if (length(precMz) == 1)
    precMz <- c(precMz - width, precMz + width)
  mz <- mz(spectrum)
  i <- intensity(spectrum)
  idx <- which(mz > precMz[1] & mz < precMz[2])
  spectrum@intensity[idx] <- 0
  return(spectrum)
}

utils.removePrecMz_list <- function(object,
                                    precMz,
                                    width = 2) {
    if (!is.numeric(precMz))
        stop("precMz must either 'NULL' or numeric.")
    if (length(precMz) > 2)
        stop("precMz must a vector of length 1 or 2.")
    if (length(precMz) == 1)
        precMz <- c(precMz - width, precMz + width)
    idx <- which(object$mz > precMz[1] & object$mz < precMz[2])
    object$int[idx] <- 0
    return(object)
}


utils.clean <- function(x, all = FALSE) {
  ## Given an numeric x, this function
  ## returns a logical b of length(x) where
  ## non-zero values and their direct
  ## 0s are TRUE so that x[b] has only
  ## non-zero values surrounded by it's
  ## original direct zero neighbours.
  ## Example
  ## x: 1 0 0 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0
  ## b: T T F T T T T T T T T T F T T T F T
  ##
  ## x[b]:  1 0 0 1 1 1 0 0 1 1 0 0 1 0 0

  n <- length(x)
  b <- rep(TRUE, n)
  if (all) {
    b[x == 0] <- FALSE
  } else {
    zeroRanges <- IRanges(sapply(x,"==",0))
    sapply(zeroRanges, function(x) {
      if (length(x) > 2)
        b[x[2:(length(x) - 1)]] <<- FALSE
    })
  }
  return(b)
}

zoom <- function(x, w = 0.05) {
    new("ReporterIons",
        mz = x,
        width = w,
        name = "xlim",
        reporterNames = paste("xlim", x, sep = "."),
        pcol = rep("grey", length(x)))
}


getBins <- function(x) {
  bins <- numeric(length(x))
  bins[1] <- 1
  for (i in 2:length(x)) {
      ifelse(x[i] == x[i-1]+1,
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
        rownames(M) <- colnames(M) <-
            paste("reporter", 1:x, sep=".")
    } else {
        if (x==4) {
            M <- matrix(c(0.929,0.059,0.002,0.000,
                          0.020,0.923,0.056,0.001,
                          0.000,0.030,0.924,0.045,
                          0.000,0.001,0.040,0.923),
                        nrow=4, byrow = TRUE)
            rownames(M) <- colnames(M) <-
                reporterNames(iTRAQ4)
        } else if (x == 6) {
            M <- matrix(c(0.939, 0.061, 0.000, 0.000, 0.000, 0.000,
                          0.005, 0.928, 0.067, 0.000, 0.000, 0.000,
                          0.000, 0.011, 0.947, 0.042, 0.000, 0.000,
                          0.000, 0.000, 0.017, 0.942, 0.041, 0.000,
                          0.000, 0.000, 0.000, 0.016, 0.963, 0.021,
                          0.000, 0.000, 0.000, 0.002, 0.032, 0.938),
                        nrow = 6, byrow = TRUE)
            rownames(M) <- colnames(M) <-
                reporterNames(TMT6)
        } else if (x == 8) {
            f <- dir(system.file("extdata", package = "MSnbase"),
                     pattern = "iTRAQ8plexPurityCorrection",
                     full.names = TRUE)
            M <- makeImpuritiesMatrix(filename = f, edit = FALSE)
            rownames(M) <- colnames(M) <- c(113:119, 121)
        } else if (x == 10) {
            ## see TMT10.R
            M <- structure(c(0.9531, 0, 0.002, 0, 0.001, 0, 0, 0, 0, 0, 0,
                             0.931, 0, 0.009, 0, 0, 0, 0, 0, 0, 0.0469, 0, 0.949, 0, 0.0053, 0, 0,
                             0, 0, 0, 0, 0.065, 0, 0.942, 0, 0.0073, 0, 0, 0, 0, 0, 0, 0.046, 0,
                             0.9678, 0, 0.013, 0, 0.001, 0, 0, 0, 0.003, 0.047, 0, 0.9678, 0,
                             0.012, 0, 0, 0, 0, 0, 0.002, 0.0259, 0, 0.962, 0, 0.029, 0, 0, 0, 0,
                             0, 0, 0.0249, 0, 0.933, 0, 0.0236, 0, 0, 0, 0, 0, 0, 0.025, 0, 0.941,
                             0, 0, 0, 0, 0, 0, 0, 0, 0.028, 0, 0.9621),
                           .Dim = c(10L, 10L),
                           .Dimnames = list(
                               c("126", "127N", "127C", "128N", "128C",
                                 "129N", "129C", "130N", "130C", "131"),
                               c("126", "127N", "127C", "128N", "128C",
                                 "129N", "129C", "130N", "130C", "131")))
        } else {
            M <- diag(x)
        }
    }
    rownames(M) <- paste("% reporter", rownames(M))
    if (edit) M <- edit(M)
    return(M)
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
  delta <- vector("list", length = length(mz.filtered))
  i <- 1
  while (length(mz.filtered) > 1) {
      m <- mz.filtered[1]
      mz.filtered <- mz.filtered[-1]
      delta[[i]] <- abs(mz.filtered-m)
      i <- i+1
  }
  return(unlist(delta))
}

utils.getMzDelta_list <- function (object, percentage) {
    idx <- order(object$int, decreasing = TRUE)
    tops <- idx[1:floor(length(idx) * percentage)]
    mz.filtered <- object$mz[tops]
    delta <- vector("list", length = length(mz.filtered))
    i <- 1
    while (length(mz.filtered) > 1) {
        m <- mz.filtered[1]
        mz.filtered <- mz.filtered[-1]
        delta[[i]] <- abs(mz.filtered - m)
        i <- i+1
    }
    return(unlist(delta))
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
  if ( is.null(dim(X)) || ncol(X) == 1 ) {
    ## vector
    unlist(mapply("[", x=split(as.vector(X), groups), i=byIdx,
                  SIMPLIFY=FALSE, USE.NAMES=FALSE))
  } else {
    ## matrix
    ans <- mapply(function(i, j) {
      X[i, , drop=FALSE][j, , drop=FALSE]
    }, i=split(1:nrow(X), groups), j=byIdx, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    do.call(rbind, ans)
  }
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

## cramer6 <- function(x, imp) {
##   if (missing(imp)) {
##     imp <- c(0, 0, 0, 6.1, 0, 0,
##              0, 0, 0.5, 6.7, 0, 0,
##              0, 0, 1.1, 4.2, 0, 0,
##              0, 0, 1.7, 4.1, 0, 0,
##              0, 0, 1.6, 2.1, 0, 0,
##              0, 0.2, 3.2, 2.8, 0, 0)
##     names(imp) <- letters[1:length(imp)]
##     impM <- matrix(imp, nrow = 6, byrow = TRUE)
##     colnames(impM) <- c("-3", "-2", "-1", "+1", "+2", "+3")
##     rownames(impM) <- 126:131
##     imp <- as.numeric(imp)
##   }
##   return(FALSE)
## }

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
##' @param f A connection object or a \code{character} string to be
##'     read in with \code{readLines(f, n = 1)}.
##' @param pattern A \code{character} string containing a regular
##'     expression to be matched to the file's header.
##' @param ... Additional parameters passed to \code{\link{strsplit}}
##'     to split the file header into individual column names.
##' @param n An \code{integer} specifying which line in file \code{f}
##'     to grep (get). Default is 1. Note that this argument must be
##'     named.
##' @return Depending on \code{value}, the matching column names of
##'     indices. In case of \code{getEcols}, a \code{character} of
##'     column names.
##' @seealso \code{\link{readMSnSet2}}
##' @author Laurent Gatto
grepEcols <- function(f, pattern, ..., n = 1)
    grep(pattern, strsplit(readLines(f, n), ...)[n][[1]])


##' @rdname grepEcols
getEcols <- function(f, ..., n = 1)
    strsplit(readLines(f, n)[n], ...)[[1]]

MSnExp.size <- function(x)
    object.size(x) + sum(unlist(unname(eapply(assayData(x),
                                              object.size))))

# convert character vector of length n to a semikolon separated character
# vector of length 1
utils.vec2ssv <- function(vec, sep=";") {
  paste0(vec, collapse=sep)
}

# convert a semikolon separated character vector of length 1 to a vector of
# length n
utils.ssv2vec <- function(ssv, sep=";", unlist=TRUE) {
  vec <- strsplit(ssv, sep)
  if (unlist) {
    return(unlist(vec))
  } else {
    return(vec)
  }
}

utils.list2ssv <- function(l, sep=";") {
  unlist(lapply(l, utils.vec2ssv, sep=sep))
}

utils.ssv2list <- function(ssv, sep=";") {
  utils.ssv2vec(ssv, unlist=FALSE, sep=sep)
}

## similar to merge(..., all.x=TRUE) but if called multiple times
## exisiting columns would not duplicated (with/without suffixes)
## but filled/overwritten using the values from y
## params: x, y, by, by.x, by.y see ?merge
##         exclude: character, columns which should excluded
##         order: logical, preserve order?
utils.leftJoin <- function(x, y, by, by.x=by, by.y=by,
                           exclude=character(), order=TRUE) {

  ## create character ids to allow ids covering several columns
  rxc <- do.call(paste, c(x[, by.x, drop=FALSE], sep=";"))
  ryc <- do.call(paste, c(y[, by.y, drop=FALSE], sep=";"))

  ## determine matching rows
  ryid <- match(rxc, ryc, 0L)
  rjid <- match(ryc, rxc, 0L)
  ryid <- ryid[ryid > 0]
  rjid <- rjid[rjid > 0]

  ## preserve order?
  if (order) {
    rjid <- sort(rjid)
  }

  cnx <- colnames(x)
  cny <- colnames(y)

  ## exclude columns
  keepx <- !cnx %in% exclude
  keepy <- !cny %in% c(exclude, by.y)
  cnx <- cnx[keepx]
  cny <- cny[keepy]
  x <- x[, keepx, drop=FALSE]
  y <- y[, keepy, drop=FALSE]

  ## start joining
  joined <- x[, cnx]

  ## only import equal columns from y
  cjid <- match(cny, cnx, 0L)
  cyid <- match(cnx, cny, 0L)

  cjid <- cjid[cjid > 0]
  cyid <- cyid[cyid > 0]

  joined[rjid, cjid] <- y[ryid, cyid]

  ## add missing columns from y
  cym <- setdiff(cny, cnx)

  if (length(cym)) {
    joined[, cym] <- NA
    joined[rjid, cym] <- y[ryid, cym]
  }

  return(joined)
}

# @param featureData fData(msexp)/fData(msset)
# @param id output of mzID::flatten(mzID(filename))
# @param fcol column name of fData data.frame used for merging
# @param icol column name of idData data.frame used for merging
# @noRd
utils.mergeSpectraAndIdentificationData <- function(featureData, id,
                                                    fcol = c("file", "acquisition.number"),
                                                    icol = c("file", "acquisitionnum")) {
  # mzR::acquisitionNum (stored in fData()[, "acquisition.number"] and
  # mzID::acquisitionnum should be identical

  if (!all(fcol %in% colnames(featureData))) {
      stop("The column(s) ", sQuote(fcol),
           " are not all in the feature data.frame!")
  }

  if (!all(icol %in% colnames(id))) {
    stop("The column(s) ", sQuote(icol),
         " are not all in the identification data.frame!")
  }

  if (sum(fcol %in% colnames(featureData)) != sum(icol %in% colnames(id))) {
    stop("The number of selected column(s) in the feature and identification ",
         "data don't match!")
  }

  ## sort id data to ensure the best matching peptide is on top in case of
  ## multiple matching peptides
  o <- do.call("order", lapply(c(icol, "rank"), function(j)id[, j]))
  id <- id[o, ]

  ## use flat version of accession/description if multiple ones are available
  id$accession <- ave(id$accession, id[, icol], FUN=utils.vec2ssv)
  id$description <- ave(id$description, id[, icol], FUN=utils.vec2ssv)

  ## remove duplicated entries
  id <- id[!duplicated(id[, icol]), ]

  featureData <- utils.leftJoin(
    x=featureData, y=id, by.x=fcol, by.y=icol,
    exclude=c("spectrumid",   # vendor specific nativeIDs
              "spectrumFile") # is stored in fileId + MSnExp@files
  )

  ## number of members in the protein group
  featureData$nprot <- sapply(utils.ssv2list(featureData$accession),
                              function(x) {
                                n <- length(x)
                                if (n == 1 && is.na(x)) return(NA)
                                n
                       })
  ## number of peptides observed for each protein
  featureData$npep.prot <- as.integer(ave(featureData$accession, featureData$pepseq, FUN=length))
  ## number of PSMs observed for each protein
  featureData$npsm.prot <- as.integer(ave(featureData$accession, featureData$accession, FUN=length))
  ## number of PSMs observed for each protein
  featureData$npsm.pep <- as.integer(ave(featureData$pepseq, featureData$pepseq, FUN=length))

  return(featureData)
}

utils.removeNoId <- function(object, fcol, keep) {
    if (!fcol %in% fvarLabels(object))
        stop(fcol, " not in fvarLabels(",
             getVariableName(match.call(), 'object'), ").")
    if (is.null(keep)) noid <- is.na(fData(object)[, fcol])
    else {
        if (!is.logical(keep))
            stop("'keep must be a logical.'")
        if (length(keep) != nrow(fData(object)))
            stop("The length of 'keep' does not match the number of spectra.")
        noid <- !keep
    }
    object <- object[!noid, ]
    object <- nologging(object, 1)
    object <- logging(object, paste0("Filtered ", sum(noid),
                                     " unidentified peptides out"))
    if (validObject(object))
        return(object)
}

utils.removeMultipleAssignment <- function(object, fcol) {
    keep <- fData(object)[, fcol] == 1
    object <- object[keep, ]
    object <- nologging(object, 1)
    object <- logging(object,
                      paste0("Removed ", sum(!keep),
                             " features assigned to multiple proteins"))
    if (validObject(object))
        return(object)
}

utils.idSummary <- function(fd) {
  if (any(!c("spectrumFile", "idFile") %in% colnames(fd))) {
    stop("No quantification/identification data found! Did you run ",
         sQuote("addIdentificationData"), "?")
  }
  idSummary <- fd[!duplicated(fd$spectrumFile), c("spectrumFile", "idFile")]
  idSummary$coverage <- sapply(idSummary$spectrumFile, function(f) {
                          round(mean(!is.na(fd$idFile[fd$spectrumFile== f])), 3)
                        })
  rownames(idSummary) <- NULL
  colnames(idSummary) <- c("spectrumFile", "idFile", "coverage")
  return(idSummary)
}

utils.removeNoIdAndMultipleAssignments <- function(object) {
    if (anyNA(fData(object)$pepseq))
        object <- removeNoId(object)
    if (any(fData(object)$nprot > 1))
        object <- removeMultipleAssignment(object)
    return(object)
}

##' Compares equality of all members of a list.
##'
##' @title Tests equality of list elements class
##' @param x A code{list}.
##' @param class A \code{character} defining the expected class.
##' @param valid A \code{logical} defining if all elements should be
##' tested for validity. Default is \code{TRUE}.
##' @return \code{TRUE} is all elements of \code{x} inherit from
##' \code{class}.
##' @author Laurent Gatto
##' @examples
##' listOf(list(), "foo")
##' listOf(list("a", "b"), "character")
##' listOf(list("a", 1), "character")
listOf <- function(x, class, valid = TRUE) {
    cla <- all(sapply(x, inherits, class))
    if (valid) val <- all(sapply(x, validObject))
    else val <- TRUE
    cla & val
}

##' Calculates a non-parametric version of the coefficient of
##' variation where the standard deviation is replaced by the median
##' absolute deviations (see \code{\link{mad}} for details) and
##' divided by the absolute value of the mean.
##'
##' Note that the \code{mad} of a single value is 0 (as opposed to
##' \code{NA} for the standard deviation, see example below).
##'
##'
##' @title Non-parametric coefficient of variation
##' @param x A \code{numeric}.
##' @param na.rm A \code{logical} (default is \code{TRUE} indicating
##' whether \code{NA} values should be stripped before the computation
##' of the median absolute deviation and mean.
##' @return A \code{numeric}.
##' @author Laurent Gatto
##' @examples
##' set.seed(1)
##' npcv(rnorm(10))
##' replicate(10, npcv(rnorm(10)))
##' npcv(1)
##' mad(1)
##' sd(1)
npcv <- function(x, na.rm = TRUE) {
    mdx <- mad(x, na.rm = na.rm)
    mdx/abs(mean(x, na.rm = na.rm))
}

##' Compares two \code{\linkS4class{MSnSet}} instances. The
##' \code{qual} and \code{processingData} slots are generally omitted.
##'
##' @title Compare two MSnSets
##' @param x First MSnSet
##' @param y Second MSnSet
##' @param qual Should the \code{qual} slots be compared? Default is
##' \code{FALSE}.
##' @param proc Should the \code{processingData} slots be compared?
##' Default is \code{FALSE}.
##' @return A \code{logical}
##' @author Laurent Gatto
compareMSnSets <- function(x, y, qual = FALSE, proc = FALSE) {
    if (!proc)  ## do not compare @processingData
        x@processingData <-
            y@processingData <- new("MSnProcess")
    if (!qual)  ## do not compare @qual
        x@qual <- y@qual
    all.equal(x, y)
}

##' Similar to colMeans but calculates the sd. Should be identical to
##' apply(x, 2, sd, na.rm).
##' based on: http://stackoverflow.com/questions/17549762/is-there-such-colsd-in-r/17551600#17551600
##' @title colSd
##' @param x matrix/data.frame
##' @param na.rm logical. Should missing values (including ‘NaN’) be omitted
##' from the calculations?
##' @return double
##' @author Sebastian Gibb <mail@@sebastiangibb.de>
##' @noRd
utils.colSd <- function(x, na.rm = TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x))
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm = na.rm) - (colMeans(x, na.rm = na.rm))^2L
  sqrt(colVar * n/(n - 1L))
}

##' Apply a function groupwise. Similar to tapply but takes a matrix as input
##' and preserve its structure and order.
##' @title applyColumnwiseByGroup
##' @param x matrix
##' @param groupBy factor/grouping index
##' @param FUN function to be applied; must work on columns, e.g. colSums
##' @param ... further arguments to FUN
##' @return modified matrix
##' @author Sebastian Gibb <mail@@sebastiangibb.de>
##' @noRd
utils.applyColumnwiseByGroup <- function(x, groupBy, FUN, ...) {
  if (!is.matrix(x)) {
    stop("x has to be a matrix!")
  }

  groupBy <- as.factor(groupBy)
  FUN <- match.fun(FUN)

  j <- split(1L:nrow(x), groupBy)
  ans <- matrix(NA_real_, nrow = nlevels(groupBy), ncol = ncol(x),
                dimnames = list(levels(groupBy), colnames(x)))

  for (i in seq(along = j)) {
    subexprs <- x[j[[i]], , drop = FALSE]
    ans[i, ] <- do.call(FUN, list(subexprs, ...))
  }

  ans
}

setMethod("trimws", "data.frame",
          function(x, which, ...) {
              for (i in 1:ncol(x)) {
                  if (inherits(x[, i], "character"))
                      x[, i] <- base::trimws(x[, i], which)
              }
              x
          })

setMethod("trimws", "MSnSet",
          function(x, which, ...) {
              fData(x) <- trimws(fData(x), which, ...)
              x <- logging(x, "Trimmed featureData white spaces")
              x
          })
