##' These functions take an instance of class
##' \code{"\linkS4class{MSnSet}"} and sets randomly selected values to
##' \code{NA}.
##'
##' \code{makeNaData} randomly selects a number \code{nNA} (or a
##' proportion \code{pNA}) of cells in the expression matrix to be set
##' to \code{NA}.
##' 
##' \code{makeNaData2} will select \code{length(nRows)} sets of rows
##' from \code{object}, each with \code{nRows[i]} rows respectively.
##' The first set will be assigned \code{nNAs[1]} missing values, the
##' second \code{nNAs[2]}, ... As opposed to \code{makeNaData}, this
##' permits to control the number of \code{NAs} per rows.
##' 
##' The \code{whichNA} can be used to extract the indices
##' of the missing values, as illustrated in the example.
##'
##' @title Create a data with missing values
##' @param object An instance of class \code{MSnSet}.
##' @param nNA The absolute number of missing values to be assigned.
##' @param pNA The proportion of missing values to be assignmed.
##' @param exclude A \code{vector} to be used to subset \code{object},
##'     defining rows that should not be used to set \code{NA}s.
##' @return An instance of class \code{MSnSet}, as \code{object}, but
##'     with the appropriate number/proportion of missing values.  The
##'     returned object has an additional feature meta-data columns,
##'     \code{nNA}
##' @author Laurent Gatto
##' @examples
##' ## Example 1
##' library(pRolocdata)
##' data(dunkley2006)
##' sum(is.na(dunkley2006))
##' dunkleyNA <- makeNaData(dunkley2006, nNA = 150)
##' processingData(dunkleyNA)
##' sum(is.na(dunkleyNA))
##' table(fData(dunkleyNA)$nNA)
##' naIdx <- whichNA(dunkleyNA)
##' head(naIdx)
##' ## Example 2
##' dunkleyNA <- makeNaData(dunkley2006, nNA = 150, exclude = 1:10)
##' processingData(dunkleyNA)
##' table(fData(dunkleyNA)$nNA[1:10])
##' table(fData(dunkleyNA)$nNA)
makeNaData <- function(object,
                       nNA, pNA,
                       exclude) {
    stopifnot(inherits(object, "MSnSet"))
    if (missing(nNA) & missing(pNA))
        stop("Provide one of 'nNA' or 'pNA'.")
    if (!missing(nNA) & !missing(pNA))
        stop("Need only one of 'nNA' or 'pNA'.")
    if (!missing(exclude)) {
        object0 <- object
        fn0 <- featureNames(object)
        if (is.logical(exclude)) {
            object <- object0[!exclude, ]
            objectX <- object0[exclude, ]
        }
        if (is.numeric(exclude)) {
            object <- object0[-(exclude), ]
            objectX <- object0[exclude, ]
        }
        if (is.character(exclude)) {
            if (!all(exclude %in% fn0))
                stop("Unknown feature names in 'exclude'")
            exl <- fn0 %in% exclude
            object <- object0[!exl, ]
            objectX <- object0[exl, ]
        }
        fData(objectX)$nNA <- 0
    }

    N <- prod(dim(object))
    if (missing(nNA)) {
        if (pNA <= 0 | pNA >= 1)
            stop("Require 0 < pNA < 1")
        nNA <- ceiling(N * pNA)
    }
    if (nNA <= 0 | nNA > N)
        stop("Require 0 < nNA > ", N)
    .nNA <- sample(N, nNA)
    .naInd <- arrayInd(.nNA, .dim = dim(object))
    .naTab <- table(.naInd[, 1])
    fData(object)$nNA <- 0
    fData(object)$nNA[as.numeric(names(.naTab))] <- .naTab
    exprs(object)[.naInd] <- NA
    stopifnot(sum(fData(object)$nNA) == sum(is.na(object)))
    msg <- {
        if (missing(pNA)) paste0("Set ", nNA, " values to NA")
        else paste0("Set ", nNA, " (", pNA, "%) values to NA")
    }
    msg <- paste(msg, date())
    if (!missing(exclude)) {
        msg <- paste0(msg, "\n  (excluding ", nrow(objectX) ," features)" )
        object <- combine(object, objectX)
        object <- object[fn0, ]
        object <- nologging(object, n = 2)
    }

    object@processingData@processing <-
        c(object@processingData@processing,
          msg)

    if (validObject(object))
        return(object)
}


##' @rdname makeNaData
##' @param nRows The number of rows for each set.
##' @param nNAs The number of missing values for each set.
##' @examples
##' ## Example 3
##' nr <- rep(10, 5)
##' na <- 1:5
##' x <- makeNaData2(dunkley2006[1:100, 1:5],
##'                  nRows = nr,
##'                  nNAs = na)
##' processingData(x)
##' (res <- table(fData(x)$nNA))
##' stopifnot(as.numeric(names(res)[-1]) ==  na)
##' stopifnot(res[-1] ==  nr)
##' ## Example 3
##' nr2 <- c(5, 12, 11, 8)
##' na2 <- c(3, 8, 1, 4)
##' x2 <- makeNaData2(dunkley2006[1:100, 1:10],
##'                   nRows = nr2,
##'                   nNAs = na2)
##' processingData(x2)
##' (res2 <- table(fData(x2)$nNA))
##' stopifnot(as.numeric(names(res2)[-1]) ==  sort(na2))
##' stopifnot(res2[-1] ==  nr2[order(na2)])
##' ## Example 5
##' nr3 <- c(5, 12, 11, 8)
##' na3 <- c(3, 8, 1, 3)
##' x3 <- makeNaData2(dunkley2006[1:100, 1:10],
##'                   nRows = nr3,
##'                   nNAs = na3)
##' processingData(x3)
##' (res3 <- table(fData(x3)$nNA))
makeNaData2 <- function(object,
                        nRows, nNAs,
                        exclude) {
  stopifnot(inherits(object, "MSnSet"))
  if (missing(nRows) | missing(nNAs))
    stop("Require 'nrows' and 'nNAs'")
  stopifnot(length(nRows) == length(nNAs))
  lNA <- length(nNAs)
  
  if (!missing(exclude)) {
    object0 <- object
    fn0 <- featureNames(object)
    if (is.logical(exclude)) {
      object <- object0[!exclude, ]
      objectX <- object0[exclude, ]
    }
    if (is.numeric(exclude)) {
      object <- object0[-(exclude), ]
      objectX <- object0[exclude, ]
    }
    if (is.character(exclude)) {
        if (!all(exclude %in% fn0))
            stop("Unknown feature names in 'exclude'")
        exl <- fn0 %in% exclude
        object <- object0[!exl, ]
        objectX <- object0[exl, ]
    }
    fData(objectX)$nNA <- 0
  }

  naRows <- sample(nrow(object), sum(nRows))
  naCols <- lapply(1:lNA, function(k) {
      replicate(nRows[k],
                sample(ncol(object), nNAs[k]))
  })
  
  for (k in 1:lNA) {
      i <- ifelse(k == 1,
                  1,
                  sum(nRows[1:(k-1)]) + 1)
      j <- sum(nRows[1:k])
      .cl <- as.numeric(naCols[[k]])
      .rw <- naRows[i:j]
      .rw <- rep(.rw, each = nNAs[k])
      stopifnot(length(.rw) == length(.cl))
      .sel <- cbind(.rw, .cl)
      exprs(object)[.sel] <- NA
  }
  
  fData(object)$nNA <-
                  colSums(apply(exprs(object), 1, is.na))
  
  msg <-  paste0("Set (", paste(nNAs, collapse = ","), ") NAs in (",
                 paste(nRows, collapse = ","), ") rows,\n  respectively ",
                 date())
  if (!missing(exclude)) {
      msg <- paste0(msg, "\n  (excluding ", nrow(objectX) ," features)" )
      object <- combine(object, objectX)
      object <- object[fn0, ]
      object <- nologging(object, n = 2)
  }

  object@processingData@processing <-
      c(object@processingData@processing, msg)

  if (validObject(object))
      return(object)
}


##' @param x A \code{matrix} or an instance of class \code{MSnSet}.
##' @rdname makeNaData
whichNA <- function(x) {
  if (inherits(x, "MSnSet"))
    x <- exprs(x)
  arrayInd(which(is.na(x)), .dim = dim(x))
}
