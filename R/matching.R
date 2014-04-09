#' similar to base::match but with tolerance
#' in contrast to match the "table" argument must be sorted and it is possible
#' to get duplicated matches
#' @param x values to be matched
#' @param table the values to be matched against, must be sorted!
#' @param nomatch the value to be returned in the case when no match is found
#' @param tolerance double, allowed deviation
#' @param relative relative (or absolute) deviation
#' @return numeric vector of the same length as "x"
relaxedMatch <- function(x, table, nomatch=NA_integer_, tolerance=25e-6,
                         relative=TRUE) {

  res <- rep(nomatch, length(x))

  if (relative) {
    if (tolerance > 1L) {
      stop(sQuote("tolerance"),
           " must be smaller than 1 for relative deviations.")
    }
    tolerance <- table*tolerance
  } else {
    tolerance <- rep_len(tolerance, length(x))
  }

  ## find left interval
  lIdx <- findInterval(x, table, rightmost.closed=FALSE, all.inside=FALSE)
  rIdx <- lIdx+1L

  ## respect borders
  lIdx[which(lIdx < 1L)] <- 1L
  rIdx[which(rIdx > length(x))] <- length(x)

  ## calculate differences for left and right
  lDiff <- abs(table[lIdx]-x)
  rDiff <- abs(table[rIdx]-x)

  potentialMatches <- ifelse(rDiff == pmin.int(lDiff, rDiff), rIdx, lIdx)
  m <- which(abs(x-table[potentialMatches]) < tolerance[potentialMatches])

  res[m] <- potentialMatches[m]
  return(res)
}

#' common peaks
#' @param x spectrum1 (MSnbase::Spectrum), to be matched
#' @param y spectrum2 (MSnbase::Spectrum), to match against
#' @param method for duplicated matches use highest/closest intensity/mz
#' @param tolerance double, allowed deviation
#' @param relative relative (or absolute) deviation
#' @return logical, TRUE for common peaks in y
commonPeaks <- function(x, y, method=c("highest", "closest"),
                         tolerance=25e-6, relative=TRUE) {
  method <- match.arg(method)

  if (peaksCount(x) == 0 || peaksCount(y) == 0) {
    return(logical(peaksCount(x)))
  }

  m <- relaxedMatch(mz(x), mz(y), nomatch=NA, tolerance=tolerance)

  if (anyDuplicated(m)) {
    if (method == "highest") {
      o <- order(intensity(x), decreasing=TRUE)
    } else {
      o <- order(abs(mz(x)-mz(y)[m]))
    }
    sortedMatches <- m[o]
    sortedMatches[which(duplicated(sortedMatches))] <- NA
    m[o] <- sortedMatches
  }

  m[which(is.na(m))] <- 0L

  return(as.logical(m))
}

#' number of common peaks
#' @param x spectrum1 (MSnbase::Spectrum)
#' @param y spectrum2 (MSnbase::Spectrum)
#' @param tolerance double, allowed deviation
#' @return double, number of common peaks
numberOfCommonPeaks <- function(x, y, tolerance=25e-6) {
  return(sum(commonPeaks(x, y, tolerance=tolerance)))
}
