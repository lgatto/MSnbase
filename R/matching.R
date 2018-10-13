#' similar to base::match but with tolerance in contrast to match the
#' 'table' argument must be sorted and it is possible to get
#' duplicated matches
#'
#' @param x values to be matched
#' @param table the values to be matched against, must be sorted!
#' @param nomatch the value to be returned in the case when no match
#' is found
#' @param tolerance double, allowed deviation
#' @param relative relative (or absolute) deviation
#' @return numeric vector of the same length as "x"
#' @noRd
relaxedMatch <- function(x, table, nomatch=NA_integer_, tolerance=25e-6,
                         relative=TRUE) {

  if (relative) {
    if (tolerance > 1L) {
      stop(sQuote("tolerance"),
           " must be smaller than 1 for relative deviations.")
    }
    tolerance <- table*tolerance
  }

  match.closest(x, table, tolerance=tolerance, nomatch=nomatch)
}

#' similar to base::match but with tolerance
#' if there are duplicated matches the highest/closest is choosen
#' @param x spectrum1 (MSnbase::Spectrum), to be matched
#' @param y spectrum2 (MSnbase::Spectrum)/double vector of mz values, to match
#' against
#' @param method for duplicated matches use "highest"/"closest" intensity/mz or
#' report "all" possible matches
#' @param tolerance double, allowed deviation
#' @param relative relative (or absolute) deviation
#' @return integer vector of the same length as "x" representing the position in
#' "y"
#' @noRd
matchPeaks <- function(x, y, method=c("highest", "closest", "all"),
                       tolerance=25e-6, relative=TRUE) {
  method <- match.arg(method)

  if (inherits(y, "Spectrum")) {
    y <- mz(y)
  }

  if (peaksCount(x) == 0 || length(y) == 0) {
    return(integer(peaksCount(x)))
  }

  m <- relaxedMatch(mz(x), y, nomatch=NA, tolerance=tolerance,
                    relative=relative)

  if (anyDuplicated(m)) {
    if (method == "highest") {
      o <- order(intensity(x), decreasing=TRUE)
    } else if (method == "closest") {
      o <- order(abs(mz(x)-y[m]))
    } else {
      o <- 1:length(x)
    }
    sortedMatches <- m[o]

    if (method != "all") {
      sortedMatches[which(duplicated(sortedMatches))] <- NA
    }
    m[o] <- sortedMatches
  }

  as.integer(m)
}

#' common peaks
#' @param x spectrum1 (MSnbase::Spectrum), to be matched
#' @param y spectrum2 (MSnbase::Spectrum), to match against
#' @param method for duplicated matches use highest/closest intensity/mz
#' @param tolerance double, allowed deviation
#' @param relative relative (or absolute) deviation
#' @return logical, TRUE for common peaks in y
#' @noRd
commonPeaks <- function(x, y, method=c("highest", "closest"),
                         tolerance=25e-6, relative=TRUE) {
  m <- matchPeaks(x, y, method=match.arg(method), tolerance=tolerance,
                  relative=relative)

  m[which(is.na(m))] <- 0L

  as.logical(m)
}

#' number of common peaks
#' @param x spectrum1 (MSnbase::Spectrum)
#' @param y spectrum2 (MSnbase::Spectrum)
#' @param tolerance double, allowed deviation
#' @param relative relative (or absolute) deviation
#' @return double, number of common peaks
#' @noRd
numberOfCommonPeaks <- function(x, y, tolerance=25e-6, relative=TRUE) {
  sum(commonPeaks(x, y, tolerance=tolerance, relative=relative))
}

#' calculate the dot product between two vectors
#'
#' Stein, S. E., and Scott, D. R. (1994).
#' Optimization and testing of mass spectral library search algorithms for
#' compound identification.
#' Journal of the American Society for Mass Spectrometry, 5(9), 859-866.
#' doi: https://doi.org/10.1016/1044-0305(94)87009-8
#'
#' Lam, H., Deutsch, E. W., Eddes, J. S., Eng, J. K., King, N., Stein, S. E.
#' and Aebersold, R. (2007)
#' Development and validation of a spectral library searching method for peptide
#' identification from MS/MS.
#' Proteomics, 7: 655-667.
#' doi: https://doi.org/10.1002/pmic.200600625
#'
#' @param x double
#' @param y double
#' @return double, length == 1
#' @noRd
dotproduct <- function(x, y) {
  as.vector(x %*% y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
}
