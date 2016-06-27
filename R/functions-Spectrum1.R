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

############################################################
##
Spectra1 <- function(peaksCount = NULL, rt=numeric(), acquisitionNum=NA_integer_,
                     scanIndex=integer(), tic=0, mz=numeric(), intensity=numeric(),
                     fromFile=integer(), centroided=FALSE, polarity=NA_integer_,
                     nvalues=integer()){
    if(length(mz) == 0 | length(intensity) == 0 | length(nvalues) == 0){
        stop("Arguments 'mz', 'intensity' and 'nvalues' are required!")
    }else{
        if(length(mz) != length(intensity))
            stop("Lengths of 'mz' and 'intensity' do not match!")
    }
    nvals <- length(nvalues)
    ## Now match all of the lengths to the length of nvalues.
    if(length(peaksCount) == 0)
        peaksCount <- nvalues
    ## rt
    if(length(rt) == 0){
        rt <- rep(NA_integer_, nvals)
    }else{
        if(length(rt) != nvals)
            stop("Length of 'rt' has to match the length of 'nvalues'!")
    }
    ## acquisitionNum
    if(length(acquisitionNum) == 1){
        rep(acquisitionNum, nvals)
    }else{
        if(length(acquisitionNum) != nvals)
            stop("Length of 'acquisitionNum' has to match the length of 'nvalues'!")
    }
    ## scanIndex
    if(length(scanIndex) == 0){
        scanIndex <- rep(NA_integer_, nvals)
    }else{
        if(length(scanIndex) != nvals)
            stop("Length of 'scanIndex' has to match the length of 'nvalues'!")
    }
    ## tic
    if(length(tic) == 1){
        rep(tic, nvals)
    }else{
        if(length(tic) != nvals)
            stop("Length of 'tic' has to match the length of 'nvalues'!")
    }
    ## fromFile
    if(length(fromFile) == 0){
        fromFile <- rep(NA_integer_, nvals)
    }else{
        if(length(fromFile) != nvals)
            stop("Length of 'fromFile' has to match the length of 'nvalues'!")
    }
    ## polarity
    if(length(polarity) == 1){
        rep(polarity, nvals)
    }else{
        if(length(polarity) != nvals)
            stop("Length of 'polarity' has to match the length of 'nvalues'!")
    }
    ## OK, now let's call C.
    res <- .Call("Multi_Spectrum1_constructor", 1L, as.integer(peaksCount),
                 as.numeric(rt),
                 as.integer(acquisitionNum),
                 as.integer(scanIndex),
                 as.numeric(tic), mz, intensity,
                 as.integer(fromFile),
                 centroided,
                 as.integer(polarity),
                 as.integer(nvalues), TRUE,
                 PACKAGE="MSnbase")
    return(res)
}



