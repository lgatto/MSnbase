show_Spectrum1 <- function(spectrum) {
  cat("Object of class \"",class(spectrum),"\"\n",sep="")
  cat(" Retention time:",formatRt(spectrum@rt),"\n")
  cat(" MSn level:",spectrum@msLevel,"\n")
  cat(" Total ion count:",spectrum@peaksCount,"\n")
  cat(" Polarity:",spectrum@polarity,"\n")
}
