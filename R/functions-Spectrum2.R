show.Spectrum2 <- function(spectrum) {
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

removeReporters.Spectrum2 <- function(object, reporters=NULL, clean=FALSE) {
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
