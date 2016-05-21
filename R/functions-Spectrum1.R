show_Spectrum1 <- function(spectrum) {
  cat("Object of class \"", class(spectrum), "\"\n", sep="")
  cat(" Retention time:", formatRt(spectrum@rt), "\n")
  cat(" MSn level:", spectrum@msLevel, "\n")
  cat(" Total ion count:", spectrum@peaksCount, "\n")
  cat(" Polarity:", spectrum@polarity, "\n")
}

############################################################
## Constructor function for Spectrum1 objects. This one uses C-code and
## is faster than a call to "new"
Spectrum1 <- function(peaksCount=length(mz), rt=numeric(),
                      acquisitionNum=NA_integer_, scanIndex=integer(), tic=0L,
                      mz=numeric(), intensity=numeric(), fromFile=integer(),
                      centroided=FALSE, polarity=NA_integer_){
    res <- .Call("Spectrum1_constructor",
                 1L, peaksCount, rt, acquisitionNum, scanIndex, tic, mz,
                 intensity, fromFile, centroided, polarity, TRUE,
                 PACKAGE="MSnbase")
    return(res)
}


