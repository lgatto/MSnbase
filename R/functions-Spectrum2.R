show.Spectrum2 <- function(spectrum) {
  cat("Object of class \"",class(spectrum),"\"\n",sep="")
  if (length(spectrum@merged)>0)
    cat(" Merged from ",length(spectrum@merged),"MSn spectra\n")
  cat(" Precursor:",spectrum@precursorMz,"\n")
  cat(" Retention time:",formatRt(spectrum@rt),"\n")
  cat(" Charge:",spectrum@precursorCharge,"\n")
  cat(" MSn level:",spectrum@msLevel,"\n")
  cat(" Peaks count:",spectrum@peaksCount,"\n")
  cat(" Total ion count:",sum(spectrum@intensity),"\n")
}
