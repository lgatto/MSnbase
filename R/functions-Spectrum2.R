#' calculate fragments from a peptide sequence for a specific Spectrum
#' only Spectrum2 is supported yet
#' @param sequence character vector of length 1
#' @param object Spectrum2 object (for Spectrum1/Spectrum object an empty
#' data.frame is returned)
#' @param tolerance double, allowed deviation for mz values to be treated as
#' equal
#' @param method matching method
#' @param relative relative (or absolute) deviation
#' @param ... further arguments passed to calculateFragments
#' @noRd
calculateFragments_Spectrum2 <- function(sequence, object, tolerance=0.1,
                                         method=c("highest", "closest", "all"),
                                         relative=FALSE, ...) {

  isValidSequence <- !missing(sequence) && !is.na(sequence) &&
                     nchar(sequence)
  isValidSpectrum <- is(object, "Spectrum2") && peaksCount(object)

  if (isValidSpectrum && isValidSequence) {

    fragments <- calculateFragments(sequence, ...)
    fragments <- fragments[order(fragments$mz), ]

    m <- matchPeaks(object, fragments$mz, tolerance=tolerance,
                    relative=relative, method=match.arg(method))
    i <- which(!is.na(m))
    fragments <- fragments[m[i], ]
    fragments$error <- fragments$mz - mz(object)[i]
    fragments$mz <- mz(object)[i]
    fragments$intensity <- intensity(object)[i]
    ## set intensity as second column
    fragments <- fragments[, c("mz", "intensity",
                               "ion", "type", "pos", "z",
                               "seq", "error")]
  } else {
    fragments <- data.frame(mz=double(), intensity=double(),
                            ion=character(), type=character(),
                            pos=integer(), z=integer(), seq=character(),
                            error=double(), stringsAsFactors=FALSE)
  }
  rownames(fragments) <- NULL
  fragments
}

show_Spectrum2 <- function(spectrum) {
  cat("Object of class \"",class(spectrum),"\"\n",sep="")
  if (length(spectrum@merged)>1)
    cat(" Merged from ",length(spectrum@merged),"MSn spectra\n")
  cat(" Precursor:",spectrum@precursorMz,"\n")
  cat(" Retention time:",formatRt(spectrum@rt),"\n")
  cat(" Charge:",spectrum@precursorCharge,"\n")
  cat(" MSn level:",spectrum@msLevel,"\n")
  cat(" Peaks count:",spectrum@peaksCount,"\n")
  cat(" Total ion count:",sum(spectrum@intensity),"\n")
}

removeReporters_Spectrum2 <- function(object, reporters=NULL, clean=FALSE) {
  ## Originally contributed by Guangchuang Yu for the plotMzDelta QC
  ## Additional modifications: setting peaks to 0 and clean argument
  ## Made removeReporters a method in version 1.1.15
  if (!is.null(reporters)) {
    mz <- mz(object)
    i <- intensity(object)
    lower <- mz(reporters) - width(reporters)
    upper <- mz(reporters) + width(reporters)
    idx <- logical(peaksCount(object))
    for (i in 1:length(lower))
      idx[mz > lower[i] & mz < upper[i]] <- TRUE
    if (sum(idx) != 0)
      object@intensity[idx] <- 0
    if (clean)
      object <- clean(object)
  }
  if (validObject(object))
    return(object)
}

