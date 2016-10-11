show_Spectrum1 <- function(spectrum) {
  cat("Object of class \"", class(spectrum), "\"\n", sep="")
  cat(" Retention time:", formatRt(spectrum@rt), "\n")
  cat(" MSn level:", spectrum@msLevel, "\n")
  cat(" Total ion count:", spectrum@peaksCount, "\n")
  cat(" Polarity:", spectrum@polarity, "\n")
}

############################################################
## Constructor function for Spectrum1 objects. This one uses C-code and
## is faster than a call to "new".
## It calls the "versioned" constructor in C that adds also the class version(s)
## (see issue #163).
Spectrum1 <- function(peaksCount = length(mz), rt = numeric(),
                      acquisitionNum = NA_integer_, scanIndex = integer(),
                      tic = 0, mz = numeric(), intensity = numeric(),
                      fromFile = integer(), centroided = NA,
                      smoothed = NA,
                      polarity = NA_integer_){
    if (tic == 0)
        tic <- sum(intensity)
    if (!is.integer(peaksCount)) peaksCount <- as.integer(peaksCount)
    if (!is.double(rt)) rt <- as.double(rt)
    if (!is.integer(acquisitionNum)) acquisitionNum <- as.integer(acquisitionNum)
    if (!is.integer(scanIndex)) scanIndex <- as.integer(scanIndex)
    if (!is.double(tic)) tic <- as.double(tic)
    if (!is.double(mz)) mz <- as.double(mz)
    if (!is.double(intensity)) intensity <- as.double(intensity)
    if (!is.integer(fromFile)) fromFile <- as.integer(fromFile)
    if (!is.logical(centroided)) centroided <- as.logical(centroided)
    if (!is.logical(smoothed)) smoothed <- as.logical(smoothed)
    if (!is.integer(polarity)) polarity <- as.integer(polarity)
    ## Define the class versions.
    versions <- list(Spectrum = getClassVersionString("Spectrum"),
                     Spectrum1 = getClassVersionString("Spectrum1"))
    res <- .Call("Spectrum1_constructor",
                 1L, peaksCount, rt, acquisitionNum, scanIndex, tic, mz,
                 intensity, fromFile, centroided, smoothed, polarity, TRUE,
                 lapply(versions, .versionToNum), PACKAGE = "MSnbase")
    return(res)
}

############################################################
## .versionToNum
## Simple helper function to convert a character version string into an
## integer vector.
.versionToNum <- function(z) {
    return(as.integer(unlist(strsplit(z, split = ".", fixed = TRUE))))
}

############################################################
## Constructor for Spectrum1 that ensures ordering of M/Z-intensity
## pairs by M/Z value (increasing).
## It calls the "versioned" constructor in C that adds also the class version(s)
## (see issue #163).
Spectrum1_mz_sorted <- function(peaksCount = length(mz), rt = numeric(),
                                acquisitionNum = NA_integer_,
                                scanIndex = integer(), tic = 0,
                                mz = numeric(), intensity = numeric(),
                                fromFile = integer(),
                                centroided = NA, smoothed = NA,
                                polarity = NA_integer_){
    if (!is.integer(peaksCount)) peaksCount <- as.integer(peaksCount)
    if (!is.double(rt)) rt <- as.double(rt)
    if (!is.integer(acquisitionNum)) acquisitionNum <- as.integer(acquisitionNum)
    if (!is.integer(scanIndex)) scanIndex <- as.integer(scanIndex)
    if (!is.double(tic)) tic <- as.double(tic)
    if (!is.double(mz)) mz <- as.double(mz)
    if (!is.double(intensity)) intensity <- as.double(intensity)
    if (!is.integer(fromFile)) fromFile <- as.integer(fromFile)
    if (!is.logical(centroided)) centroided <- as.logical(centroided)
    if (!is.logical(smoothed)) smoothed <- as.logical(smoothed)
    if (!is.integer(polarity)) polarity <- as.integer(polarity)
    ## Define the class versions.
    versions <- list(Spectrum = getClassVersionString("Spectrum"),
                     Spectrum1 = getClassVersionString("Spectrum1"))
    res <- .Call("Spectrum1_constructor_mz_sorted",
                 1L, peaksCount, rt, acquisitionNum, scanIndex, tic, mz,
                 intensity, fromFile, centroided, smoothed, polarity, TRUE,
                 lapply(versions, .versionToNum), PACKAGE="MSnbase")
    return(res)
}


############################################################
## Constructor to create multiple Spectrum1 objects ensuring that
## M/Z-intensity pairs are ordered by M/Z.
## It calls the "versioned" constructor in C that adds also the class version(s)
## (see issue #163).
Spectra1_mz_sorted <- function(peaksCount = NULL, rt = numeric(),
                               acquisitionNum = NA_integer_,
                               scanIndex = integer(), tic = 0, mz = numeric(),
                               intensity = numeric(), fromFile = integer(),
                               centroided = NA, smoothed = NA,
                               polarity = NA_integer_, nvalues = integer()) {
    if (length(mz) == 0 | length(intensity) == 0 | length(nvalues) == 0) {
        stop("Arguments 'mz', 'intensity' and 'nvalues' are required!")
    } else {
        if (length(mz) != length(intensity))
            stop("Lengths of 'mz' and 'intensity' do not match!")
    }
    nvals <- length(nvalues)
    ## Now match all of the lengths to the length of nvalues.
    if (length(peaksCount) == 0)
        peaksCount <- nvalues
    ## rt
    if (length(rt) == 0){
        rt <- rep(NA_integer_, nvals)
    } else {
        if (length(rt) != nvals)
            stop("Length of 'rt' has to match the length of 'nvalues'!")
    }
    ## acquisitionNum
    if (length(acquisitionNum) == 1) {
        acquisitionNum <- rep(acquisitionNum, nvals)
    } else {
        if (length(acquisitionNum) != nvals)
            stop("Length of 'acquisitionNum' has to match the length of 'nvalues'!")
    }
    ## scanIndex
    if (length(scanIndex) == 0){
        scanIndex <- rep(NA_integer_, nvals)
    } else {
        if (length(scanIndex) != nvals)
            stop("Length of 'scanIndex' has to match the length of 'nvalues'!")
    }
    ## tic
    if (length(tic) == 1){
        tic <- rep(tic, nvals)
    } else {
        if (length(tic) != nvals)
            stop("Length of 'tic' has to match the length of 'nvalues'!")
    }
    ## fromFile
    if (length(fromFile) == 0){
        fromFile <- rep(NA_integer_, nvals)
    } else {
        if (length(fromFile) != nvals)
            stop("Length of 'fromFile' has to match the length of 'nvalues'!")
    }
    ## polarity
    if (length(polarity) == 1){
        polarity <- rep(polarity, nvals)
    } else {
        if (length(polarity) != nvals)
            stop("Length of 'polarity' has to match the length of 'nvalues'!")
    }
    ## centroided and smoothed
    if (length(centroided) == 1){
        centroided <- rep(centroided, nvals)
    } else {
        if (length(centroided) != nvals)
            stop("Length of 'centroided' has to match length of 'nvalues'!")
    }
    if (length(smoothed) == 1){
        smoothed <- rep(smoothed, nvals)
    } else {
        if (length(smoothed) != nvals)
            stop("Length of 'smoothed' has to match length of 'nvalues'!")
    }
    ## Ensure the arguments are in correct format
    if (!is.integer(peaksCount)) peaksCount <- as.integer(peaksCount)
    if (!is.double(rt)) rt <- as.double(rt)
    if (!is.integer(acquisitionNum)) acquisitionNum <- as.integer(acquisitionNum)
    if (!is.integer(scanIndex)) scanIndex <- as.integer(scanIndex)
    if (!is.double(tic)) tic <- as.double(tic)
    if (!is.double(mz)) mz <- as.double(mz)
    if (!is.double(intensity)) intensity <- as.double(intensity)
    if (!is.integer(fromFile)) fromFile <- as.integer(fromFile)
    if (!is.logical(centroided)) centroided <- as.logical(centroided)
    if (!is.logical(smoothed)) smoothed <- as.logical(smoothed)
    if (!is.integer(polarity)) polarity <- as.integer(polarity)
    ## Define the class versions.
    versions <- list(Spectrum = getClassVersionString("Spectrum"),
                     Spectrum1 = getClassVersionString("Spectrum1"))
    ## OK, now let's call C.
    res <- .Call("Multi_Spectrum1_constructor_mz_sorted",
                 1L,
                 peaksCount,
                 rt,
                 acquisitionNum,
                 scanIndex,
                 tic, mz, intensity,
                 fromFile,
                 centroided,
                 smoothed,
                 polarity,
                 as.integer(nvalues), TRUE,
                 lapply(versions, .versionToNum),
                 PACKAGE = "MSnbase")
    return(res)
}
