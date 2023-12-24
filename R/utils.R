##'  This function is used to convert retention times. Conversion is
##'  seconds to/from the more human friendly format "mm:sec". The
##'  implementation is from [MsCoreUtils::formatRt()].
##'
##' @title Format Retention Time
##'
##' @param rt retention time in seconds (`numeric`) or "mm:sec"
##'     (`character`).
##'
##' @return A vector of same length as `rt`.
##'
##' @author Laurent Gatto and Sebastian Gibb
##'
##' @examples
##'
##' formatRt(1524)
##' formatRt("25:24")
formatRt <- function(rt)
    MsCoreUtils::formatRt(rt)


#' @param fileIds numeric, could be a vector
#' @param spectrumIds numeric, could be a vector
#' @param nFiles numeric, max number of files
#' @param nSpectra numeric, max number of spectra
#' @return character, in the format F001.S0001
#' @noRd
formatFileSpectrumNames <- function(fileIds, spectrumIds,
                                    nFiles=length(fileIds),
                                    nSpectra=length(spectrumIds)) {
  digits <- ceiling(log10(c(nFiles, nSpectra) + 1L))

  if (length(fileIds) != 1L && length(spectrumIds) != length(fileIds)) {
    stop("Length of 'fileIds' has to be one or equal to ",
         "the length of 'spectrumIds'.")
  }

  sprintf(paste0("F%0", digits[1L], "d.S%0", digits[2L], "d"),
          fileIds, spectrumIds)
}

utils.removePeaks_centroided <- function(int, t) {
    int[int <= t] <- 0L
    int
}

utils.removePeaks <- function(int, t) {
    peakRanges <- as(int > 0L, "IRanges")
    toLow <- max(extractList(int, peakRanges)) <= t
    replaceROWS(int, peakRanges[toLow], 0L)
}

## For internal use - use utils.removePrecMz_Spectrum that will set
## the paramters based on data accessed directly in the spectrum
## object.
utils.removePrecMz <- function(mz, int, precMz, tolerance = 25e-6) {
    if (!is.numeric(precMz) || length(precMz) != 1L) {
        stop("precMz must be numeric of length 1.")
    }

    i <- relaxedMatch(precMz, mz, tolerance = tolerance)

    if (!is.na(i)) {
        peakRanges <- as(int > 0L, "IRanges")
        i <- findOverlaps(IRanges(i, width = 1L), peakRanges,
                          type = "within", select = "first")
        if (!is.na(i)) {
            int <- replaceROWS(int, peakRanges[i], 0L)
        }
    }
    int
}

utils.removePrecMz_Spectrum <- function(spectrum,
                                        precMz = NULL,
                                        tolerance = 25e-6) {
    if (is.null(precMz))
        precMz <- precursorMz(spectrum)
    if (!is.numeric(precMz))
        stop("precMz must either 'NULL' or numeric.")
    spectrum@intensity <- utils.removePrecMz(mz(spectrum),
                                             intensity(spectrum),
                                             precMz = precMz,
                                             tolerance = tolerance)
    spectrum
}

utils.removePrecMz_list <- function(object, precMz, tolerance = 25e-6) {
    idx <- which(object$mz > precMz[1] & object$mz < precMz[2])
    object$int <- utils.removePrecMz(object$mz,
                                     object$int,
                                     precMz = precMz,
                                     tolerance = tolerance)
    object
}

#' Removes zeros from input except the ones that in the direct neighbourhood of
#' non-zero values.
#'
#' @param x \code{numeric}, vector to be cleaned
#' @param all \code{logical}, should all zeros be removed?
#' @param na.rm \code{logical}, should NAs removed before looking for zeros?
#' @return logical vector, \code{TRUE} for keeping the value
#' @note The return value for \code{NA} is always \code{FALSE}.
#' @examples
#' x <- c(1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0)
#' #      T, T, F, T, T, T, T, T, T, T, T, F, F
#' r <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
#'        FALSE, FALSE)
#' stopifnot(utils.clean(x) == r)
#' @noRd
utils.clean <- function(x, all=FALSE, na.rm=FALSE) {
  notNA <- !is.na(x)
  notZero <- x != 0 & notNA

  if (all) {
    notZero
  } else if (na.rm) {
    notNA[notNA] <- utils.enableNeighbours(notZero[notNA])
    notNA
  } else {
    utils.enableNeighbours(notZero)
  }
}

#' Switch FALSE to TRUE in the direct neighborhod of TRUE.
#' (used in utils.clean)
#'
#' @param x logical
#' @return logical
#' @examples
#' x <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
#'        FALSE, FALSE)
#' r <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
#'        FALSE, FALSE)
#' stopifnot(utils.enableNeighbours(x) == r)
#' @noRd
utils.enableNeighbours <- function(x) {
  stopifnot(is.logical(x))
  x | c(x[-1], FALSE) | c(FALSE, x[-length(x)])
}

zoom <- function(x, w = 0.05) {
    new("ReporterIons",
        mz = x,
        width = w,
        name = "xlim",
        reporterNames = paste("xlim", x, sep = "."),
        pcol = rep("grey", length(x)))
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
        rownames(M) <- colnames(M) <- rownames(m)
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
            M <- structure(c(0.95, 0, 0.003, 0, 0, 0, 0, 0, 0, 0, 0, 0.94, 0,
                             0.004, 0, 0, 0, 0, 0, 0, 0.05, 0, 0.949, 0, 0.006,
                             0, 0, 0, 0, 0, 0, 0.058, 0, 0.955, 0, 0.008, 0,
                             0.001, 0, 0, 0, 0, 0.048, 0, 0.964, 0, 0.014, 0, 0,
                             0, 0, 0, 0, 0.041, 0, 0.957, 0, 0.015, 0, 0.002, 0,
                             0, 0, 0, 0.03, 0, 0.962, 0, 0.017, 0, 0, 0, 0, 0,
                             0, 0.035, 0, 0.928, 0, 0.02, 0, 0, 0, 0, 0, 0,
                             0.024, 0, 0.965, 0, 0, 0, 0, 0, 0, 0, 0, 0.024, 0,
                             0.956),
                           .Dim = c(10L, 10L),
                           .Dimnames = list(
                               c("126", "127N", "127C", "128N", "128C", "129N",
                                 "129C", "130N", "130C", "131"),
                               c("126", "127N", "127C", "128N", "128C", "129N",
                                 "129C", "130N", "130C", "131")))
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

#' rowwise max, similar to rowwise mean via rowMeans
#'
#' @param x matrix
#' @param na.rm logical
#' @return double vector with maximum values per row
#' @seealso Biobase::rowMax (could not handle missing values/NA)
#' @noRd
.rowMaxs <- function(x, na.rm=FALSE) {
  stopifnot(is.matrix(x))
  if (na.rm) {
    x[is.na(x)] <- -Inf
  }
  nr <- nrow(x)
  x[(max.col(x, ties.method="first") - 1L) * nr + 1L:nr]
}

#' summarise rows by an user-given function
#'
#' @param x matrix
#' @param fun function to summarise rows, if \code{fun} equals
#' \code{sum}/\code{mean} the more efficient \code{rowSums}/\code{rowMeans} are
#' used.
#' @param ... further arguments passed to \code{fun}
#' @return double, summarised rows
#' @noRd
.summariseRows <- function(x, fun, ...) {
  stopifnot(is.matrix(x))
  stopifnot(is.function(fun))

  if (identical(fun, sum)) {
    rowSums(x, ...)
  } else if (identical(fun, mean)) {
    rowMeans(x, ...)
  } else {
    apply(x, 1L, fun, ...)
  }
}

#' find top n indices of each group
#'
#' @param x matrix
#' @param groupBy factor/character of length \code{nrow(x)}
#' @param n consider just the top \code{n} values
#' @param fun function to summarise rows
#' @param ... further arguments passed to \code{fun}
#' @return double, indices sorted by summarising function \code{fun}
#' @noRd
.topIdx <- function(x, groupBy, n, fun, ...) {
  if (n < 1) {
    stop(sQuote("n"), " has to be greater or equal than 1.")
  }
  if (nrow(x) != length(groupBy)) {
    stop(sQuote("nrow(x)"), " and ", sQuote("length(groupBy)"),
         " have to be equal.")
  }
  rs <- .summariseRows(x, fun, ...)
  o <- order(as.double(rs), decreasing=TRUE, na.last=TRUE)
  idx <- unlist(lapply(split(o, groupBy[o]), "[", 1:n), use.names=FALSE)
  idx[!is.na(idx)]
}

## Computes header from assay data by-passing cache
.header <- function(object) {
  if (length(object) == 0)
    return(data.frame())
  if (all(msLevel(object) == 1)) {
    ln <- length(object)
    nas <- rep(NA, ln)
    hd <- list(fileIdx = fromFile(object),
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
    hd <- list(fileIdx = fromFile(object),
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
  if (missing(pattern)) {
    stop(sQuote("pattern"), " must not be missing.")
  }
  if (!is.matrix(x)) {
    stop(sQuote("x"), " must be a matrix.")
  }
  if (nchar(pattern) != ncol(x)) {
    stop("The ", sQuote("pattern"), " must be equal to the number of columns.")
  }
  pattern <- strsplit(pattern, "")[[1L]]
  if (!all(unique(pattern) %in% c("0", "1"))) {
    stop(sQuote("pattern"), " must be composed of '0' or '1' defining columns",
         " with or without 'NA's.")
  }
  pattern == "1"
}

getRowsFromPattern <- function(x, pattern) {
  cols <- getColsFromPattern(x, pattern)
  x <- x[, cols, drop=FALSE]
  rowSums(is.na(x)) == 0
}

.filterNA <- function(x, pNA=0L, pattern) {
  if (!is.matrix(x)) {
    stop(sQuote("x"), " must be a matrix.")
  }
  if (!is.numeric(pNA)) {
    stop(sQuote("pNA"), " must be numeric.")
  }
  if (length(pNA) > 1) {
    stop(sQuote("pNA"), " must be of length one.")
  }

  if (missing(pattern)) { ## using pNA
    if (pNA > 1) {
      pNA <- 1
    }
    if (pNA < 0) {
      pNA <- 0
    }
    k <- rowSums(is.na(x)) / ncol(x)
    k <= pNA
  } else { ## using pattern
    getRowsFromPattern(x, pattern)
  }
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
    msg <- paste0(msg, " [", date(), "]")
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

## convert vector of length n to a semicolon separated character
## vector of length 1
utils.vec2ssv <- function(vec, sep=";") {
    paste0(vec, collapse=sep)
}

## converts an n by m data.frame into an 1 by m data.frame where the
## vector columns of length n are converted to a semicolon separated
## character vector of length 1
utils.vec2ssv.data.frame <- function(x, sep = ";", exclude) {
    if (nrow(x) == 1L)
        return(x)
    nms <- names(x)
    if (missing(exclude)) {
        ans <- lapply(x, utils.vec2ssv)
    } else {
        if (is.numeric(exclude)) {
            x0 <- x[, exclude, drop = FALSE]
            x <- x[, -exclude, drop = FALSE]
        } else if (is.character(exclude)) {
            ## convert to logical, making sure that the column names
            ## to be excluded are present
            stopifnot(all(exclude %in% names(x)))
            exclude <- names(x) %in% exclude
            x0 <- x[, exclude, drop = FALSE]
            x <- x[, !exclude, drop = FALSE]
        } else if (is.logical(exclude)) {
            x0 <- x[, exclude, drop = FALSE]
            x <- x[, !exclude, drop = FALSE]
        } else {
            stop("Can only exclude numeric, characters or logicals.")
        }
        ans <- lapply(x, utils.vec2ssv)
        x0 <- lapply(x0, head, n = 1L)
        ans <- c(x0, ans)
        ans <- ans[nms] ## preserve original order
    }
    data.frame(ans, stringsAsFactors = FALSE)
}

## convert a semicolon separated character vector of length 1 to a
## vector of length n
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

## @param featureData fData(msexp)/fData(msset)
## @param id output of mzID::flatten(mzID(filename))
## @param fcol column name of fData data.frame used for merging
## @param icol column name of idData data.frame used for merging
## @noRd
utils.mergeSpectraAndIdentificationData <- function(featureData, id,
                                                    fcol, icol, acc,
                                                    desc, pepseq,
                                                    rank = "rank") {
    ## mzR::acquisitionNum (stored in fData()[, "acquisition.number"] and
    ## mzID::acquisitionnum should be identical
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
    o <- do.call("order", lapply(c(icol, rank), function(j) id[, j]))
    id <- id[o, ]

    ## use flat version of accession/description if multiple ones are available
    id[, acc] <- ave(as.character(id[, acc]), id[, icol], FUN = utils.vec2ssv)
    id[, desc] <- ave(as.character(id[, desc]), id[, icol], FUN = utils.vec2ssv)

    ## remove duplicated entries
    id <- id[!duplicated(id[, icol]), ]

    featureData <- utils.leftJoin(
        x = featureData, y = id, by.x = fcol, by.y = icol,
        exclude = c("spectrumid",   # vendor specific nativeIDs
                    "spectrumID",
                    "spectrumFile") # is stored in fileId + MSnExp@files
    )

    ## number of members in the protein group
    featureData$nprot <- sapply(utils.ssv2list(featureData[, acc]),
                                function(x) {
                                    n <- length(x)
                                    if (n == 1 && is.na(x)) return(NA)
                                    n
                                })

    ## number of peptides observed for each protein
    featureData$npep.prot <- as.integer(ave(featureData[, acc],
                                            featureData[, pepseq],
                                            FUN = length))

    ## number of PSMs observed for each protein
    featureData$npsm.prot <- as.integer(ave(featureData[, acc],
                                            featureData[, acc],
                                            FUN = length))

    ## number of PSMs observed for each protein
    featureData$npsm.pep <- as.integer(ave(featureData[, pepseq],
                                           featureData[, pepseq],
                                           FUN = length))

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

utils.removeMultipleAssignment <- function(object, nprot = "nprot") {
    keep <- which(fData(object)[, nprot] == 1)
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
        round(mean(!is.na(fd$idFile[fd$spectrumFile == f])), 3)
    })
    rownames(idSummary) <- NULL
    colnames(idSummary) <- c("spectrumFile", "idFile", "coverage")
    return(idSummary)
}

utils.removeNoIdAndMultipleAssignments <-
    function(object, pepseq = "sequence", nprot = "nprot") {
        if (anyNA(fData(object)[, pepseq]))
            object <- removeNoId(object, pepseq)
        if (any(fData(object)[, nprot] > 1))
            object <- removeMultipleAssignment(object, nprot)
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

##' Similar to rowsum but calculates the mean. It is slower than colMeans but
##' supports grouping variables. See ?rowsum for details.
##' @param x matrix
##' @param group a vector/factor of grouping
##' @param reorder if TRUE the rows are ordered by `sort(unique(group))`
##' @param na.rm logical. Should missing values (including `NaN`) be omitted
##' @return matrix
##' @author Sebastian Gibb <mail@@sebastiangibb.de>
##' @noRd
rowmean <- function(x, group, reorder=FALSE, na.rm=FALSE) {
  if (na.rm) {
    nna <- !is.na(x)
    mode(nna) <- "numeric"
  } else {
    nna <- x
    nna[] <- 1
  }
  nna <- rowsum(nna, group=group, reorder=reorder, na.rm=na.rm)
  rs <- rowsum(x, group=group, reorder=reorder, na.rm=na.rm)
  rs/nna
}

##' Similar to rowsum but calculates the sd.
##' See ?rowsum for details.
##' @param x matrix
##' @param group a vector/factor of grouping
##' @param reorder if TRUE the rows are ordered by `sort(unique(group))`
##' @param na.rm logical. Should missing values (including `NaN`) be omitted
##' @return matrix
##' @author Sebastian Gibb <mail@@sebastiangibb.de>
##' @noRd
rowsd <- function(x, group, reorder=FALSE, na.rm=FALSE) {
  if (na.rm) {
    nna <- !is.na(x)
    mode(nna) <- "numeric"
  } else {
    nna <- x
    nna[] <- 1
  }
  nna <- rowsum(nna, group=group, reorder=reorder, na.rm=na.rm)
  nna[nna == 1] <- NA_real_            # return NA if n == 1 (similar to sd)
  var <- rowmean(x*x, group=group, reorder=reorder, na.rm=na.rm) -
    rowmean(x, group=group, reorder=reorder, na.rm=na.rm)^2L
  sqrt(var * nna/(nna - 1L))
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

setMethod("isEmpty", "environment",
          function(x) length(ls(x)) == 0)

## Simple helper to help differentiate between on disk and in
## memory objects.
isOnDisk <- function(object)
    any(grepl("spectraProcessingQueue", slotNames(object)))

## Simple function to determine whether parallel or serial processing should be
## performed.
## Check testthat/test_OnDiskMSnExp_benchmarks.R for performance comparisons.
## Parameter object is expected to beb a
getBpParam <- function(object, BPPARAM=bpparam()) {
    parallel_thresh <- options()$MSnbase$PARALLEL_THRESH
    if (is.null(parallel_thresh) )
        parallel_thresh <- 1000
    ## If it's empty, return SerialParam
    if (length(object) == 0)
        return(SerialParam())
    if (length(object) < parallel_thresh)
        return(SerialParam())
    return(BPPARAM)
}

countAndPrint <- function(x) {
    if (length(x) == 0)
        return("")
    tb <- table(x)
    paste(paste0(names(tb), " (", tb, ")"), collapse = ", ")
}

## see issue #131
.isCentroided <- function(pk, k = 0.025, qtl = 0.9) {
        .qtl <- quantile(pk[, 2], qtl)
        x <- pk[pk[, 2] > .qtl, 1]
        quantile(diff(x), 0.25) > k
}

##' @title Reads profile/centroided mode from an mzML file
##' @param x An instance of \code{MSnExp} or \code{OnDiskMSnExp}
##' @return A \code{logical}
##' @noRd
.isCentroidedFromFile <- function(f) {
    if (!requireNamespace("XML"))
        stop("Please install the XML package to use this functionality.")
    xml <- XML::xmlParse(f)
    x <- XML::xpathSApply(xml,
                          "//x:spectrum/x:cvParam[@accession='MS:1000127' or @accession='MS:1000128']/../@index |
                    //x:cvParam[@accession='MS:1000127' or @accession='MS:1000128']/@name",
                    namespaces = c(x = "http://psi.hupo.org/ms/mzml"))
    index <- as.double(x[seq(1, length(x), by = 2)])
    res <- rep(NA, length(index))
    res[grepl("centroid", x[seq(2, length(x), by = 2)])] <- TRUE
    res[grepl("profile",  x[seq(2, length(x), by = 2)])] <- FALSE
    res
}

## Returns the extension of the file. If that extension is on of the
## usual archive extensions, as defined in gexts, then the last part
## after the dot is removed and the extension is extracted again.
.fileExt <- function(f,
                     gexts = c("gz", "gzip", "bzip", "bzip2", "xz",
                               "zip")) {
    ext <- tools::file_ext(f)
    if (ext %in% gexts) {
        f <- basename(f)
        f <- sub("\\.[a-z]+$", "", f)
        ext <- .fileExt(f)
    }
    ext
}

## see this comment
## https://github.com/lgatto/MSnbase/issues/183#issuecomment-273512931
## for some background about this function
.firstMsLevel <- function(object) {
  if (inherits(object, "OnDiskMSnExp")) msLevel(object)[1]
  else msLevel(object[[1]])
}

#' @title Open an MS file using the mzR package
#'
#' @description Opens an MS file using the mzR package determining the corrent
#'     backend based on the file ending of the specified file.
#'
#' @param x \code{character(1)}: the file name.
#'
#' @return A file handle to the opened MS file.
#'
#' @author Johannes Rainer
#'
#' @noRd
.openMSfile <- function(x) {
    if (missing(x) || length(x) != 1)
        stop("parameter 'x' has to be of length 1")
    mzR::openMSfile(x, backend = NULL)
}


##' This function produces the opposite as the \code{stringsAsFactors}
##' argument in the \code{data.frame} or \code{read.table} functions;
##' it converts \code{factors} columns to \code{characters}.
##'
##' @title Converts factors to strings
##' @param x A \code{data.frame}
##' @return A \code{data.frame} where \code{factors} are converted to
##'     \code{characters}.
##' @author Laurent Gatto
##' @examples
##' data(iris)
##' str(iris)
##' str(factorsAsStrings(iris))
factorsAsStrings <- function(x) {
    x <- lapply(x,
                   function(xx) {
                       if (is.factor(xx)) as.character(xx)
                       else xx
                   })
    data.frame(x, stringsAsFactors = FALSE)
}

##' Convert a \code{vector} of characters to camel case by replacing
##' dots by captial letters.
##'
##' @title Convert to camel case by replacing dots by captial letters
##' @param x A \code{vector} to be transformed to camel case.
##' @param prefix An optional \code{character} of length one. Any
##'     additional elements are ignores.
##' @return A \code{character} of same length as \code{x}.
##' @author Laurent Gatto
##' @examples
##' nms <- c("aa.foo", "ab.bar")
##' makeCamelCase(nms)
##' makeCamelCase(nms, prefix = "x")
makeCamelCase <- function(x, prefix) {
    if (!missing(prefix))
        x <- paste(prefix[1], x, sep = ".")
    gsub('\\.(\\w?)', '\\U\\1', x, perl = TRUE)
}


##' Reduce a data.frame so that the (primary) key column contains only
##' unique entries and other columns pertaining to that entry are
##' combined into semicolon-separated values into a single
##' row/observation.
##'
##' An important side-effect of reducing a `data.frame` is that all
##' columns other than the key are converted to characters when they
##' are collapsed to a semi-column separated value (even if only one
##' value is present) as soon as one observation of transformed.
##'
##' @title Reduce a data.frame
##' @param x A \code{data.frame}.
##' @param key The column name (currenly only one is supported) to be
##'     used as primary key.
##' @param sep The separator. Default is \code{;}.
##' @return A reduced \code{data.frame}.
##' @author Laurent Gatto
##' @examples
##' dfr <- data.frame(A = c(1, 1, 2),
##'                   B = c("x", "x", "z"),
##'                   C = LETTERS[1:3])
##' dfr
##' dfr2 <- reduce(dfr, key = "A")
##' dfr2
##' ## column A used as key is still num
##' str(dfr2)
##' dfr3 <- reduce(dfr, key = "B")
##' dfr3
##' ## A is converted to chr; B remains factor
##' str(dfr3)
##' dfr4 <- data.frame(A = 1:3,
##'                    B = LETTERS[1:3],
##'                    C = c(TRUE, FALSE, NA))
##' ## No effect of reducing, column classes are maintained
##' str(reduce(dfr4, key = "B"))
setMethod("reduce", "data.frame",
          function(x, key, sep = ";") {
              if (nrow(x) %in% c(0, 1))
                  return(x)
              if (missing(key))
                  stop("Need a key column to reduce the data.frame")
              if (length(key) != 1L)
                  stop("Key must be of length 1.")
              if (!key %in% names(x))
                  stop("key not found in column names.")
              ans <- by(x, x[, key], utils.vec2ssv.data.frame, exclude = key)
              ans <- do.call(rbind, ans)
              rownames(ans) <- NULL
              ans
          })

.reduce_list <- function(x) {
    x <- x[lengths(x) > 0]
    sel <- sapply(x, function(xx) any(xx != ""))
    x[sel]
}

#' @param fd data.frame, feature data (columns required: acquisitionNum,
#' precursorScanNum)
#' @param an integer, acquisitionNum of spectrum of interest (parent and
#' children will be selected)
#' @noRd
.filterSpectraHierarchy <- function(fd, an) {
    if (!is.data.frame(fd)) {
        stop("'fd' is not a data.frame")
    }
    if (!all(c("acquisitionNum", "precursorScanNum") %in% colnames(fd))) {
        stop("column(s) acquisitionNum/precursorScanNum is/are missing")
    }
    ## we could use recursion which is slow in R
    ## or reformat the adjacency list into a nested tree
    ## list model but most ms data are limited to at most 3 levels and the
    ## filtering isn't done very often, so we use for loops here

    parents <- logical(nrow(fd))

    ## find current scan
    parents[fd$acquisitionNum %in% an] <- TRUE
    children <- parents

    ## find parent scan
    nLastParents <- 0L
    nParents <- 1L
    while (nLastParents < nParents) {
        parents[fd$acquisitionNum %in% fd$precursorScanNum[parents]] <- TRUE
        nLastParents <- nParents
        nParents <- sum(parents)
    }

    ## find children scans
    nLastChildren <- 0L
    nChildren <- 1L
    while (nLastChildren < nChildren) {
        children[fd$precursorScanNum %in% fd$acquisitionNum[children]] <- TRUE
        nLastChildren <- nChildren
        nChildren <- sum(children)
    }
    parents | children
}

windowIndices <- function(i, hws, n) {
    stopifnot(i <= n)
    max(1L, i - hws):min(n, i + hws)
}

#' The function aggregates `x` for `toBin` falling into bins defined
#' by `breaks` using the `fun` function.
#'
#' @details
#'
#' This is a combination of the code from the former bin_Spectrum.
#'
#' @param x `numeric` with the values that should be binned.
#'
#' @param toBin `numeric`, same length than `x`, with values to be used for the
#'     binning.
#'
#' @param binSize `numeric(1)` with the size of the bins.
#'
#' @param breaks `numeric` defining the breaks/bins.
#'
#' @param fun `function` to be used to aggregate values of `x` falling into the
#'     bins defined by `breaks`.
#'
#' @return `list` with elements `x` and `mids` being the aggregated values
#'     of `x` for values in `toBin` falling within each bin and the bin mid
#'     points.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @noRd
.bin_values <- function(x, toBin, binSize = 1, breaks = seq(floor(min(toBin)),
                                                            ceiling(max(toBin)),
                                                            by = binSize),
                        fun = max) {
    if (length(x) != length(toBin))
        stop("lengths of 'x' and 'toBin' have to match.")
    fun <- match.fun(fun)
    breaks <- .fix_breaks(breaks, range(toBin))
    nbrks <- length(breaks)
    idx <- findInterval(toBin, breaks)
    ## Ensure that indices are within breaks.
    idx[which(idx < 1L)] <- 1L
    idx[which(idx >= nbrks)] <- nbrks - 1L

    ints <- double(nbrks - 1L)
    ints[unique(idx)] <- unlist(lapply(base::split(x, idx), fun),
                                use.names = FALSE)
    list(x = ints, mids = (breaks[-nbrks] + breaks[-1L]) / 2L)
}

#' Simple function to ensure that breaks (for binning) are span al leat the
#' expected range.
#'
#' @param brks `numeric` with *breaks* such as calculated by `seq`.
#'
#' @param rng `numeric(2)` with the range of original numeric values on which
#'     the breaks were calculated.
#'
#' @noRd
.fix_breaks <- function(brks, rng) {
    ## Assuming breaks being sorted.
    if (brks[length(brks)] <= rng[2])
        brks <- c(brks, max((rng[2] + 1e-6),
                            brks[length(brks)] + mean(diff(brks))))
    brks
}

##' Helper functions to check whether raw files contain spectra or
##' chromatograms.
##'
##' @title Checks if raw data files have any spectra or chromatograms
##' @param files A `character()` with raw data filenames.
##' @return A `logical(n)` where `n == length(x)` with `TRUE` if that
##'     files contains at least one spectrum, `FALSE` otherwise.
##' @author Laurent Gatto
##' @rdname hasSpectraOrChromatograms
##' @md
##' @examples
##' f <- msdata::proteomics(full.names = TRUE)[1:2]
##' hasSpectra(f)
##' hasChromatograms(f)
hasSpectra <- function(files) {
    sapply(files, mzR:::.hasSpectra)
}

##' @rdname hasSpectraOrChromatograms
hasChromatograms <- function(files) {
    sapply(files, mzR:::.hasChromatograms)
}

#' @title Get the index of the particular element for each level
#'
#' `levelIndex` returns the index of the first, middle or last element for
#' each level of a factor within the factor.
#'
#' @param x `factor` or `vector` that can be converted into a `factor`
#'
#' @param which `character` defining for which element the index should be
#'     returned, can be either `"first"`, `"middle"` or `"last"`.
#'
#' @return `integer` same length than `levels(x)` with the index for each
#'     level in `x`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' f <- factor(c("a", "a", "b", "a", "b", "c", "c", "b", "d", "d", "d"))
#' f
#'
#' levelIndex(f, which = "first")
#' levelIndex(f, which = "middle")
#' levelIndex(f, which = "last")
#'
#' f <- factor(c("a", "a", "b", "a", "b", "c", "c", "b", "d", "d", "d"),
#'     levels = c("d", "a", "c", "b"))
#' levelIndex(f, which = "first")
#' levelIndex(f, which = "middle")
#' levelIndex(f, which = "last")
levelIndex <- function(x, which = c("first", "middle", "last")) {
    x <- as.factor(x)
    res <- switch(match.arg(which),
                  "first" = match(levels(x), x),
                  "last" = length(x) - match(levels(x), rev(x)) + 1L,
                  "middle" = vapply(levels(x), function(z) {
                      idx <- which(x == z)
                      idx[ceiling(length(idx) / 2L)]
                  }, integer(1), USE.NAMES = FALSE))
    names(res) <- levels(x)
    res
}
