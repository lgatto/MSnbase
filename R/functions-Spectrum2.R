#' calculate fragments from a peptide sequence for a specific Spectrum
#' only Spectrum2 is supported yet
#' @param sequence character vector of length 1
#' @param object Spectrum2 object (for Spectrum1/Spectrum object an empty
#' data.frame is returned)
#' @param tolerance double, allowed deviation for mz values to be treated as
#' equal
#' @param method matching method
#' @param relative relative (or absolute) deviation
#' @param ... further arguments passed to `PSMarch:::.calculateFragments()`
#' @noRd
calculateFragments_Spectrum2 <- function(sequence, object, tolerance=0.1,
                                         method=c("highest", "closest", "all"),
                                         relative=FALSE, ...) {

  isValidSequence <- !missing(sequence) && !is.na(sequence) &&
                     nchar(sequence)
  isValidSpectrum <- is(object, "Spectrum2") && peaksCount(object)

  if (isValidSpectrum && isValidSequence) {

    fragments <- calculateFragments(sequence, ...)
    fragments <- fragments[base::order(fragments$mz), ]

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
    if (length(spectrum@merged) > 1)
        cat(" Merged from ",length(spectrum@merged),"MSn spectra\n")
    cat(" Precursor:",spectrum@precursorMz,"\n")
    if (length(spectrum@rt))
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


############################################################
## C-level constructor for multiple Spectrum2 instances (i.e. returns a
## list of instances).
## This enforces ordering of M/Z-intensity value pairs by M/Z (in C).
## It also uses the "versioned" constructor in C that adds also the class
## version(s) (see issue #163).
Spectra2_mz_sorted <- function(peaksCount = NULL, rt = numeric(),
                               acquisitionNum = integer(),
                               scanIndex = integer(), tic = numeric(),
                               mz = numeric(),
                               intensity = numeric(), fromFile = integer(),
                               centroided = logical(), smoothed = logical(),
                               polarity = integer(), msLevel = 2L,
                               merged = numeric(), precScanNum = integer(),
                               precursorMz = numeric(),
                               precursorIntensity = numeric(),
                               precursorCharge = integer(),
                               collisionEnergy = numeric(),
                               nvalues = integer()) {
    ## Argument check; make sure all is OK before calling C.
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
    ## Initialize empty vectors.
    emptyInt <- rep(NA_integer_, nvals)
    emptyNum <- as.numeric(emptyInt)
    emptyLog <- as.logical(emptyInt)
    ## Now match all of the lengths to the length of nvalues.
    if (!length(peaksCount))
        peaksCount <- nvalues
    ## rt
    if (!length(rt)) {
        rt <- emptyNum
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
    ## centroided
    if (!length(centroided)) {
        centroided <- emptyLog
    } else {
        if (length(centroided) != nvals)
            stop("Length of 'centroided' has to match the length of 'nvalues'!")
    }
    ## smoothed
    if (!length(smoothed)) {
        smoothed <- emptyLog
    } else {
        if (length(smoothed) != nvals)
            stop("Length of 'smoothed' has to match the length of 'nvalues'!")
    }
    ## msLevel
    if (!length(msLevel)) {
        msLevel <- emptyInt
    } else {
        if (length(msLevel) != nvals)
            stop("Length of 'msLevel' has to match the length of 'nvalues'!")
    }
    ## merged
    if (!length(merged)) {
        merged <- emptyNum
    } else {
        if (length(merged) != nvals)
            stop("Length of 'merged' has to match the length of 'nvalues'!")
    }
    ## precScanNum
    if (!length(precScanNum)) {
        precScanNum <- emptyInt
    } else {
        if (length(precScanNum) != nvals)
            stop("Length of 'precScanNum' has to match the length of 'nvalues'!")
    }
    ## precursorMz
    if (!length(precursorMz)) {
        precursorMz <- emptyNum
    } else {
        if (length(precursorMz) != nvals)
            stop("Length of 'precursorMz' has to match the length of 'nvalues'!")
    }
    ## precursorIntensity
    if (!length(precursorIntensity)) {
        precursorIntensity <- emptyNum
    } else {
        if (length(precursorIntensity) != nvals)
            stop("Length of 'precursorIntensity' has to match the length of 'nvalues'!")
    }
    ## precursorCharge
    if (!length(precursorCharge)) {
        precursorCharge <- emptyNum
    } else {
        if (length(precursorCharge) != nvals)
            stop("Length of 'precursorCharge' has to match the length of 'nvalues'!")
    }
    ## collisionEnergy
    if (!length(collisionEnergy)) {
        collisionEnergy <- emptyNum
    } else {
        if (length(collisionEnergy) != nvals)
            stop("Length of 'collisionEnergy' has to match the length of 'nvalues'!")
    }
    ## Ensure that we have the correct data types before passing to C
    if (!is.integer(msLevel)) msLevel <- as.integer(msLevel)
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
    if (!is.double(merged)) merged <- as.double(merged)
    if (!is.integer(precScanNum)) precScanNum <- as.integer(precScanNum)
    if (!is.double(precursorMz)) precursorMz <- as.double(precursorMz)
    if (!is.double(precursorIntensity)) precursorIntensity <- as.double(precursorIntensity)
    if (!is.integer(precursorCharge)) precursorCharge <- as.integer(precursorCharge)
    if (!is.double(collisionEnergy)) collisionEnergy <- as.double(collisionEnergy)
    ## Define the class versions.
    versions <- list(Spectrum = getClassVersionString("Spectrum"),
                     Spectrum2 = getClassVersionString("Spectrum2"))
    ## OK, now let's call C.
    res <- .Call("Multi_Spectrum2_constructor_mz_sorted",
                 msLevel,
                 peaksCount,
                 rt,
                 acquisitionNum,
                 scanIndex,
                 tic, mz, intensity,
                 fromFile,
                 centroided,
                 smoothed,
                 polarity,
                 merged,
                 precScanNum,
                 precursorMz,
                 precursorIntensity,
                 precursorCharge,
                 collisionEnergy,
                 as.integer(nvalues), TRUE,
                 lapply(versions, .versionToNum),
                 PACKAGE = "MSnbase")
    return(res)
}

############################################################
## Constructor function for Spectrum2 objects. This one uses C-code and
## is faster than a call to "new"
## It calls the "versioned" constructor in C that adds also the class version(s)
## (see issue #163).
Spectrum2 <- function(msLevel = 2L, peaksCount = length(mz), rt = numeric(),
                      acquisitionNum = NA_integer_, scanIndex = integer(),
                      tic = 0L, mz = numeric(), intensity = numeric(),
                      fromFile = integer(), centroided = NA,
                      smoothed = NA, polarity = NA_integer_,
                      merged = 1, precScanNum = NA_integer_,
                      precursorMz = NA, precursorIntensity = NA,
                      precursorCharge = NA_integer_, collisionEnergy = NA) {
    if (tic == 0)
        tic <- sum(intensity)
    ## Ensure that we have the correct data types before passing to C
    if (!is.integer(msLevel)) msLevel <- as.integer(msLevel)
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
    if (!is.double(merged)) merged <- as.double(merged)
    if (!is.integer(precScanNum)) precScanNum <- as.integer(precScanNum)
    if (!is.double(precursorMz)) precursorMz <- as.double(precursorMz)
    if (!is.double(precursorIntensity)) precursorIntensity <- as.double(precursorIntensity)
    if (!is.integer(precursorCharge)) precursorCharge <- as.integer(precursorCharge)
    if (!is.double(collisionEnergy)) collisionEnergy <- as.double(collisionEnergy)
    ## Define the class versions.
    versions <- list(Spectrum = getClassVersionString("Spectrum"),
                     Spectrum2 = getClassVersionString("Spectrum2"))
    res <- .Call("Spectrum2_constructor",
                 msLevel, peaksCount, rt, acquisitionNum, scanIndex, tic, mz,
                 intensity, fromFile, centroided, smoothed, polarity,
                 merged, precScanNum, precursorMz, precursorIntensity,
                 precursorCharge, collisionEnergy,
                 TRUE, lapply(versions, .versionToNum),
                 PACKAGE = "MSnbase")
    return(res)
}

############################################################
## Constructor function for Spectrum2 objects. This one uses C-code and
## is faster than a call to "new"
## It calls the "versioned" constructor in C that adds also the class version(s)
## (see issue #163).
Spectrum2_mz_sorted <- function(msLevel = 2L, peaksCount = length(mz), rt = numeric(),
                                acquisitionNum = NA_integer_, scanIndex = integer(),
                                tic = 0L, mz = numeric(), intensity = numeric(),
                                fromFile = integer(), centroided = NA,
                                smoothed = NA, polarity = NA_integer_,
                                merged = 1, precScanNum = NA_integer_,
                                precursorMz = NA, precursorIntensity = NA,
                                precursorCharge = NA_integer_, collisionEnergy = NA) {
    ## if (tic == 0)
    ##     tic <- sum(intensity)
    ## Ensure that we have the correct data types before passing to C
    if (!is.integer(msLevel)) msLevel <- as.integer(msLevel)
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
    if (!is.double(merged)) merged <- as.double(merged)
    if (!is.integer(precScanNum)) precScanNum <- as.integer(precScanNum)
    if (!is.double(precursorMz)) precursorMz <- as.double(precursorMz)
    if (!is.double(precursorIntensity)) precursorIntensity <- as.double(precursorIntensity)
    if (!is.integer(precursorCharge)) precursorCharge <- as.integer(precursorCharge)
    if (!is.double(collisionEnergy)) collisionEnergy <- as.double(collisionEnergy)
    ## Define the class versions.
    versions <- list(Spectrum = getClassVersionString("Spectrum"),
                     Spectrum2 = getClassVersionString("Spectrum2"))
    .Call("Spectrum2_constructor_mz_sorted",
          msLevel, peaksCount, rt, acquisitionNum, scanIndex, tic, mz,
          intensity, fromFile, centroided, smoothed, polarity,
          merged, precScanNum, precursorMz, precursorIntensity,
          precursorCharge, collisionEnergy,
          TRUE, lapply(versions, .versionToNum),
          PACKAGE = "MSnbase")
}
