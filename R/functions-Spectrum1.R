show_Spectrum1 <- function(spectrum) {
    cat("Object of class \"", class(spectrum), "\"\n", sep="")
    if (length(spectrum@rt))
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
    .Call("Spectrum1_constructor_mz_sorted",
          1L, peaksCount, rt, acquisitionNum, scanIndex, tic, mz,
          intensity, fromFile, centroided, smoothed, polarity, TRUE,
          lapply(versions, .versionToNum), PACKAGE="MSnbase")
}


############################################################
## Constructor to create multiple Spectrum1 objects ensuring that
## M/Z-intensity pairs are ordered by M/Z.
## It calls the "versioned" constructor in C that adds also the class version(s)
## (see issue #163).
Spectra1_mz_sorted <- function(peaksCount = NULL, rt = numeric(),
                               acquisitionNum = integer(),
                               scanIndex = integer(), tic = numeric(),
                               mz = numeric(),
                               intensity = numeric(), fromFile = integer(),
                               centroided = logical(), smoothed = logical(),
                               polarity = integer(), nvalues = integer()) {
    ## Fix issue #215: remove check to allow empty spectra
    ## if (length(mz) == 0 | length(intensity) == 0 | length(nvalues) == 0) {
    ##     stop("Arguments 'mz', 'intensity' and 'nvalues' are required!")
    ## }
    if (sum(nvalues) != length(mz))
        stop("Length of 'mz' does not match with the number of values per ",
             "spectrum")
    if (length(mz) != length(intensity))
        stop("Lengths of 'mz' and 'intensity' do not match!")
    nvals <- length(nvalues)
    emptyInt <- rep(NA_integer_, nvals)
    emptyNum <- as.numeric(emptyInt)
    emptyLog <- as.logical(emptyInt)
    ## Now match all of the lengths to the length of nvalues.
    if (!length(peaksCount))
        peaksCount <- nvalues
    ## rt
    if (!length(rt)) {
        rt <- emptyInt
    } else {
        if (length(rt) != nvals)
            stop("Length of 'rt' has to match the length of 'nvalues'!")
    }
    ## acquisitionNum
    if (!length(acquisitionNum)) {
        acquisitionNum <- emptyInt
    } else {
        if (length(acquisitionNum) != nvals)
            stop("Length of 'acquisitionNum' has to match the length of 'nvalues'!")
    }
    ## scanIndex
    if (!length(scanIndex)) {
        scanIndex <- emptyInt
    } else {
        if (length(scanIndex) != nvals)
            stop("Length of 'scanIndex' has to match the length of 'nvalues'!")
    }
    ## tic
    if (!length(tic)) {
        tic <- emptyNum
    } else {
        if (length(tic) != nvals)
            stop("Length of 'tic' has to match the length of 'nvalues'!")
    }
    ## fromFile
    if (!length(fromFile)) {
        fromFile <- emptyInt
    } else {
        if (length(fromFile) != nvals)
            stop("Length of 'fromFile' has to match the length of 'nvalues'!")
    }
    ## polarity
    if (!length(polarity)) {
        polarity <- emptyInt
    } else {
        if (length(polarity) != nvals)
            stop("Length of 'polarity' has to match the length of 'nvalues'!")
    }
    ## centroided and smoothed
    if (!length(centroided)) {
        centroided <- emptyLog
    } else {
        if (length(centroided) != nvals)
            stop("Length of 'centroided' has to match length of 'nvalues'!")
    }
    if (!length(smoothed)) {
        smoothed <- emptyLog
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
    .Call("Multi_Spectrum1_constructor_mz_sorted",
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
}
