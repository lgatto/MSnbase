############################################################
## Methods for OnDiskMSnExp objects.

############################################################
## initialize
##
setMethod("initialize",
          signature(.Object = "OnDiskMSnExp"),
          function(.Object, ...){
              callNextMethod()
          })

setMethod("header",
          c("OnDiskMSnExp", "missing"),
          function(object) {
              .Deprecated("fData")
              fData(object)
          })

setMethod("header",
          c("OnDiskMSnExp", "numeric"),
          function(object, scans) {
              .Deprecated("fData")
              fData(object)[scans, ]
          })


setReplaceMethod("featureNames",
                 c("OnDiskMSnExp", "ANY"),
                 function(object, value) {
                     fd <- featureData(object)
                     featureNames(fd) <- value
                     object@featureData <- fd
                     if (validObject(object))
                         return(object)
                 })


as.MSnExp.OnDiskMSnExp <- function(x, ...)
    as(x, "OnDiskMSnExp")

setAs("OnDiskMSnExp", "MSnExp",
      function(from) {
          if (length(unique(msLevel(from))) != 1) {
              msg <- c("'MSnExp' objects don't support multiple MS levels. ",
                       "Please use 'filterMsLevel' to keep only one level ",
                       "before coercing. ")
              stop(paste(msg, collapse = "\n"))
          }
          ans <- new("MSnExp")
          slts <- intersect(slotNames("OnDiskMSnExp"),
                            slotNames("MSnExp"))
          slts <- setdiff(slts,
                          c(".__classVersion__", ".cache", "assayData"))
          for (sl in slts)
              slot(ans, sl) <- slot(from, sl)
          ans@assayData <- list2env(spectra(from))
          ## Remove column fileIdx from the fData (issue #308).
          fData(ans) <- fData(ans)[, colnames(fData(ans)) != "fileIdx"]
          lockEnvironment(ans@assayData, bindings = TRUE)
          ans <- logging(ans, "Coerced from OnDiskMSnExp")
          if (validObject(ans))
              return(ans)
      })



############################################################
## processingQueue
processingQueue <- function(object) {
    return(object@spectraProcessingQueue)
}

############################################################
## msLevel
##
## Extract the msLevel info for all spectra in an OnDiskMSnExp
## object. In contrast to the MSnExp we're not getting that
## from the individual Spectrum objects in @assayData, but
## from the featureData.
setMethod("msLevel", "OnDiskMSnExp", function(object) {
    msl <- featureData(object)$msLevel
    ## That should not happen!
    if (is.null(msl))
        stop("The 'OnDiskMSnExp' does not have msLevel information in the featureData!")
    names(msl) <- featureNames(object)
    return(msl)
})

############################################################
## fromFile
##
## Extract the index from which file the spectra in the data set
## derive from.
setMethod("fromFile", "OnDiskMSnExp", function(object) {
    fidx <- fData(object)$fileIdx
    names(fidx) <- featureNames(object)
    return(fidx)
})
setReplaceMethod("fromFile", signature(object = "OnDiskMSnExp",
                                       value = "integer"),
                 function(object, value) {
                     if (length(object) != length(value))
                         stop("Length of replacement value is different from the number of spectra.")
                     object@featureData$fileIdx <- value
                     valMsg <- validObject(object)
                     if (valMsg) {
                         return(object)
                     } else {
                         stop(valMsg)
                     }
                 })

############################################################
## length
##
## Get the length, i.e. the number of spectra we've got (from the
## featureData).
setMethod("length", "OnDiskMSnExp", function(x) {
    return(nrow(fData(x)))
})

############################################################
## scanIndex
##
## Get the scan index for each spectrum in each file. We're extracting
## that from the featureData.
setMethod("scanIndex", "OnDiskMSnExp", function(object) {
    scIdx <- fData(object)$spIdx
    names(scIdx) <- featureNames(object)
    return(scIdx)
})

############################################################
## precScanNum
##
## Get the precursor scan index for each spectrum in each file. We're
## extracting that from the featureData.
setMethod("precScanNum", "OnDiskMSnExp",
          function(object) {
              pscan <- fData(object)$precursorScanNum
              names(pscan) <- featureNames(object)
              return(pscan)
          })

############################################################
## precursorIntensity
##
## Get the precursor intensity for each spectrum in each file. We're
## extracting that from the featureData.
setMethod("precursorIntensity", "OnDiskMSnExp",
          function(object) {
              pint <- fData(object)$precursorIntensity
              names(pint) <- featureNames(object)
              return(pint)
          })

############################################################
## acquisitionNum
##
## Get the acquisition number for each spectrum in each file. We're extracting
## that from the featureData.
setMethod("acquisitionNum", "OnDiskMSnExp",
          function(object) {
              aIdx <- fData(object)$acquisitionNum
              names(aIdx) <- featureNames(object)
              return(aIdx)
          })

############################################################
## centroided
##
## Getter/setter for the centroided information; extracting this from
## the featureData.
setMethod("centroided", "OnDiskMSnExp",
          function(object, na.fail = FALSE) {
              val <- fData(object)$centroided
              if (na.fail & any(is.na(val)))
                  stop("Mode is undefined. See ?isCentroided for details.",
                       call. = FALSE)
              names(val) <- featureNames(object)
              return(val)
          })

setReplaceMethod("centroided",
                 signature(object = "OnDiskMSnExp", value = "logical"),
                 function(object, msLevel., ..., value) {
                     if (missing(msLevel.)) {
                         if (length(value) == 1)
                             value <- rep(value, length(object))
                         if (length(object) != length(value))
                             stop("Length of replacement value is different than number of spectra.")
                         fData(object)$centroided <- value
                     } else {
                         sel <- fData(object)$msLevel == msLevel.
                         fData(object)$centroided[sel] <- value
                     }
                     if (validObject(object))
                         return(object)
                 })

setMethod("isCentroided", "OnDiskMSnExp",
          function(object, ..., verbose = isMSnbaseVerbose()) {
              res <- spectrapply(object, function(z, ...)
                  .isCentroided(as(z, "data.frame"), ...))
              ctrd <- unlist(res, use.names = FALSE)
              if (verbose) print(table(ctrd, msLevel(object)))
              names(ctrd) <- featureNames(object)
              ctrd
          })

############################################################
## smoothed
##
## Getter/setter for the centroided information; extracting this from
## the featureData.
setMethod("smoothed", "OnDiskMSnExp",
          function(object){
              val <- fData(object)$smoothed
              names(val) <- featureNames(object)
              return(val)
          })

setReplaceMethod("smoothed",
                 signature(object = "OnDiskMSnExp", value = "logical"),
                 function(object, msLevel., ..., value) {
                     if (missing(msLevel.)) {
                         if (length(value) == 1)
                             value <- rep(value, length(object))
                         if (length(object) != length(value))
                             stop("Length of replacement value is different than number of spectra.")
                         fData(object)$smoothed <- value
                     } else {
                         sel <- fData(object)$msLevel == msLevel.
                         fData(object)$smoothed[sel] <- value
                     }
                     if (validObject(object))
                         return(object)
                 })

############################################################
## rtime
##
## Get the retention time
setMethod("rtime", "OnDiskMSnExp",
          function(object){
              vals <- fData(object)$retentionTime
              names(vals) <- featureNames(object)
              return(vals)
          })

############################################################
## tic
##
## Get the total ion Current (from the featureData, thus it will not
## be recalculated).
setMethod("tic", "OnDiskMSnExp",
          function(object, initial = TRUE, BPPARAM = bpparam()) {
              if (initial) {
                  vals <- fData(object)$totIonCurrent
              } else {
                  ## Calculate the value.
                  vals <- unlist(
                      spectrapply(object, FUN = function(z) sum(z@intensity),
                                  BPPARAM = BPPARAM))
              }
              names(vals) <- featureNames(object)
              vals
          })

############################################################
## bpi
##
## Get the base peak intensity from each spectrum. If initial = TRUE
## the basePeakIntensity from the feaureData is returned, i.e. the value
## extracted from the header of the raw mzML file. If initial = FALSE the
## actual base peak intensity will be calculated for each spectrum by first
## applying all eventual data manipulation operations.
setMethod("bpi", "OnDiskMSnExp",
          function(object, initial = TRUE, BPPARAM = bpparam()) {
              if (initial) {
                  vals <- fData(object)$basePeakIntensity
              } else {
                  ## Calculate the value.
                  vals <- unlist(
                      spectrapply(object, FUN = function(z) max(z@intensity),
                                  BPPARAM = BPPARAM))
              }
              names(vals) <- featureNames(object)
              vals
          })


############################################################
## collisionEnergy
##
## Get the collision Energy (from the featureData, thus it will not be
## recalculated).
setMethod("collisionEnergy", "OnDiskMSnExp",
          function(object) {
              vals <- fData(object)$collisionEnergy
              names(vals) <- featureNames(object)
              return(vals)
          })

############################################################
## polarity
##
## Get the polarity.
setMethod("polarity", "OnDiskMSnExp",
          function(object){
              vals <- fData(object)$polarity
              names(vals) <- featureNames(object)
              return(vals)
          })

#############################################################
## peaksCount
##
## Extract the peaksCount data; if the spectraProcessingQueue is empty
## we're just returning the peaksCount from the featureData, otherwise we
## check the function in there and eventually load the data and apply it.
setMethod("peaksCount",
          signature(object = "OnDiskMSnExp", scans = "missing"),
          function(object, scans, BPPARAM = bpparam()) {
              ## The feature data contains the original peaks
              ## count. This method fetches the peaks count from the
              ## (possibly processed) spectra.
              ##
              ## Add methods here that would not require raw data reading.
              skipFun <- c("removePeaks")

              if (length(object@spectraProcessingQueue) > 0) {
                  recalc <- any(unlist(lapply(object@spectraProcessingQueue,
                                              function(z) {
                                                  if (any(z@FUN == skipFun))
                                                      return(FALSE)
                                                  return(TRUE)
                                              })))
              } else {
                  ## No need to calculate the peaks count; we can use the
                  ## information from the feature data.
                  recalc <- FALSE
              }
              ## An important point here is that we DON'T want to get
              ## all of the data from all files in one go; that would
              ## require eventually lots of memory! It's better to do
              ## that per file; that way we could also do that in
              ## parallel.
              if (recalc) {
                  vals <- unlist(spectrapply(object,
                                             FUN = peaksCount,
                                             BPPARAM = BPPARAM))
              } else {
                  vals <- fData(object)$originalPeaksCount
              }
              names(vals) <- featureNames(object)
              return(vals)
          })

############################################################
## ionCount
##
## Calculate the ion count, i.e. the sum of intensities per spectrum.
setMethod("ionCount", "OnDiskMSnExp",
          function(object, BPPARAM = bpparam()) {
              ## An important point here is that we DON'T want to get
              ## all of the data from all files in one go; that would
              ## require eventually lots of memory! It's better to do
              ## that per file; that way we could also do that in
              ## parallel.
              vals <- spectrapply(object,
                                  FUN = function(y) return(sum(y@intensity)),
                                  BPPARAM = BPPARAM)
              return(unlist(vals))
          })

############################################################
## spectra
##
## Extract the spectra of an experiment by reading the raw data from
## the original files, applying processing steps from the queue.
setMethod("spectra",
          "OnDiskMSnExp",
          function(object, BPPARAM = bpparam()){
              return(spectrapply(object, BPPARAM = BPPARAM))
          })

############################################################
## assayData
##
## Read the full data, put it into an environment and return that.
setMethod("assayData", "OnDiskMSnExp",
          function (object) {
              return(list2env(spectra(object)))
          })

############################################################
## intensity
##
## Extract the intensity values of individual spectra. This means we
## have to read all of the data, create Spectrum objects, apply
## eventual processing steps and return the intensities.
setMethod("intensity", "OnDiskMSnExp",
          function(object, BPPARAM = bpparam())
              return(spectrapply(object, FUN = intensity, BPPARAM = BPPARAM)))


############################################################
## mz
##
## Extract the mz values of individual spectra.
setMethod("mz", "OnDiskMSnExp",
          function(object, BPPARAM = bpparam())
    return(spectrapply(object, FUN = mz, BPPARAM = BPPARAM)))

############################################################
## [[
##
## Extract individual spectra; allow extraction of multiple
## spectra too, i.e. length i > 1.
setMethod("[[", "OnDiskMSnExp",
          function(x, i, j = "missing", drop = "missing") {
              if (length(i) != 1)
                  stop("subscript out of bounds")
              if (is.character(i))
                  i <- base::match(i, featureNames(x))
              if (any(is.na(i)))
                  stop("subscript out of bounds")
              ## Keep the fromFile information; that would be
              ## changed by [ subsetting.
              fromF <- fromFile(x)[[i]]
              ## Subsetting first the OnDiskMSnExp with subsequent
              ## spectra extraction is considerably faster.
              spctr <- spectra(x[i])[[1]]
              ## Re-setting the original fromFile values.
              spctr@fromFile <- fromF
              return(spctr)
          })

setMethod("quantify",
          signature = signature("OnDiskMSnExp"),
          function(object,
                   method = c(
                       "trapezoidation", "max", "sum",
                       "SI", "SIgi", "SIn",
                       "SAF", "NSAF",
                       "count"),
                   reporters,
                   wd = width(reporters),
                   strict = FALSE,
                   BPPARAM,
                   verbose = isMSnbaseVerbose(),
                   ...) {
              method <- match.arg(method)
              if (!all(msLevel(object) >= 2)) {
                  message("Currently only MS > 1 quantitation: filtering MS2 spectra\n.",
                          "Use filterMsLevel() to filter appropriate levels.")
                  object <- filterMsLevel(object, msLevel. = 2L)
                  if (length(object) == 0L)
                      stop("Empty MSnExp data.")
              }
              ## MS2 isobaric
              if (method %in% c("trapezoidation", "max", "sum")) {
                  if (!inherits(reporters, "ReporterIons"))
                      stop("Argument 'reporters' must inherit from ",
                           "'ReporterIons' class.")
                  if (missing(BPPARAM)) {
                      BPPARAM <- bpparam()
                      if (verbose)
                          message("Using default parallel backend: ",
                                  class(BPPARAM)[1])
                  }
                  if (method != "max")
                      stop("Not yet implemented - see issue #130")
                  if (!verbose)
                      suppressMessages(quantify_OnDiskMSnExp_max(object,
                                                                 reporters,
                                                                 wd,
                                                                 BPPARAM))
                  else quantify_OnDiskMSnExp_max(object, reporters,
                                                 wd, BPPARAM)
              } else if (method == "count") {
                  count_MSnSet(object)
              } else {
                  ## the following assumes that the appropriate fcols
                  ## are available
                  object <- utils.removeNoIdAndMultipleAssignments(object)
                  if (method %in% c("SI", "SIgi", "SIn")) SI(object, method, ...)
                  else SAF(object, method, ...)
              }
          })


############################################################
## [
##
## Subset by [
setMethod("[", "OnDiskMSnExp",
          function(x, i, j, ..., drop = TRUE) {
              if (!missing(j))
                  stop("Subsetting by column ('j =", j, "') is not supported")
              if (!(is.logical(i) | is.numeric(i)))
                  stop("'i' has to be numeric or logical")
              if (is.logical(i))
                  i <- which(i)
              i <- base::sort(i)  ## Force sorted!
              ## Now subset the featureData. The function will
              ## complain if i is outside the range.
              if (length(i) == 0) x@featureData <- x@featureData[0, ]
              else x@featureData <- subsetFeatureDataBy(featureData(x),
                                                        index = i)
              ## Check also that processingData and experimentData
              ## match the *new* featureData:
              file <- base::sort(unique(fromFile(x)))
              ## o Sub-set the phenoData.
              pd <- phenoData(x)[file, , drop = FALSE]
              pData(pd) <- droplevels(pData(pd))
              x@phenoData <- pd
              ## Sub-set the files.
              x@processingData@files <- x@processingData@files[file]
              ## Sub-set the experimentData:
              ## o instrumentManufacturer
              ## o instrumentModel
              ## o ionSource
              ## o analyser
              ## o detectorType
              expD <- experimentData(x)
              if (length(expD@instrumentManufacturer))
                  expD@instrumentManufacturer <- expD@instrumentManufacturer[file]
              else expD@instrumentManufacturer <- expD@instrumentManufacturer
              if (length(expD@instrumentModel))
                  expD@instrumentModel <- expD@instrumentModel[file]
              else expD@instrumentModel <- expD@instrumentModel
              if (length(expD@ionSource))
                  expD@ionSource <- expD@ionSource[file]
              else expD@ionSource <- expD@ionSource
              if (length(expD@analyser))
                  expD@analyser <- expD@analyser[file]
              else expD@analyser <- expD@ionSource
              if (length(expD@detectorType))
                  expD@detectorType <- expD@detectorType[file]
              else expD@detectorType <- expD@ionSource
              x@experimentData <- expD
              ## Update fromFile in spectra/featureData.
              fromFile(x) <- base::match(fromFile(x), file)
              if (validObject(x))
                  return(x)
          })

############################################################
## precursorMz
setMethod("precursorMz", "OnDiskMSnExp",
          function(object) {
              return(precursorValue_OnDiskMSnExp(object,
                                                 column = "precursorMZ"))
          })

############################################################
## precursorCharge
setMethod("precursorCharge", "OnDiskMSnExp",
          function(object) {
              return(precursorValue_OnDiskMSnExp(object,
                                                 column = "precursorCharge"))
          })

############################################################
## precursorIntensity
setMethod("precursorIntensity", "OnDiskMSnExp",
          function(object) {
              return(precursorValue_OnDiskMSnExp(object,
                                                 column = "precursorIntensity"))
          })

############################################################
## precScanNum
setMethod("precScanNum", "OnDiskMSnExp",
          function(object) {
              return(precursorValue_OnDiskMSnExp(object,
                                                 column = "precursorScanNum"))
          })

############################################################
## extractPrecSpectra
##
## Just subset the featureData based on the prec.
setMethod("extractPrecSpectra", signature = signature(object = "OnDiskMSnExp",
                                                      prec = "numeric"),
          function(object, prec){
              ## Match prec with the precursorMz.
              gotEm <- which(precursorMz(object) %in% prec)
              if (length(gotEm) == 0)
                  stop("None of the precursor MZ values found in",
                       " 'OnDiskMSnExp' object.")
              if (length(gotEm) != length(prec))
                  warning(length(prec) - length(gotEm), " precursor MZ",
                          " values not found in 'OnDiskMSnExp' object.")
              ## Subset the object; sorting the indices!
              object <- object[base::sort(gotEm)]
              ## Add processingData.
              object@processingData@processing <-
                  c(object@processingData@processing,
                    paste(length(prec), " (", length(gotEm),
                          ") precursors (spectra)", " extracted: ",
                          date(), sep = ""))
              return(object)
          })

##============================================================
##  --  DATA MANIPULATION METHODS
##
##------------------------------------------------------------

############################################################
## removePeaks
##
## Add a "removePeaks" ProcessingStep to the queue and update
## the processingData information of the object.
setMethod("removePeaks", signature("OnDiskMSnExp"),
          function(object, t = "min", verbose = isMSnbaseVerbose(), msLevel.) {
              if (missing(t))
                  t <- "min"
              if (!is.numeric(t)) {
                  if (t != "min")
                      stop("Argument 't' has to be either numeric or 'min'!")
              }
              if (missing(msLevel.)) {
                  msLevel. <- base::sort(unique(msLevel(object)))
              } else {
                  if (!is.numeric(msLevel.))
                      stop("'msLevel' must be numeric!")
              }
              ps <- ProcessingStep("removePeaks", list(t = t,
                                                       msLevel. = msLevel.))
              ## Append the processing step to the queue.
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@removedPeaks <-
                  c(object@processingData@removedPeaks,
                    as.character(t))
              object@processingData@processing <-
                  c(object@processingData@processing,
                    paste0("Curves <= ", t, " in MS level(s) ",
                           paste0(msLevel., collapse = ", "),
                           "set to '0': ", date()))
              return(object)
          })

############################################################
## clean
##
## Add a "clean" ProcessingStep to the queue and update
## the processingData information of the object.
setMethod("clean", signature("OnDiskMSnExp"),
          function(object, all = FALSE, verbose = isMSnbaseVerbose(), msLevel.) {
              if (!is.logical(all))
                  stop("Argument 'all' is supposed to be a logical!")
              if (missing(msLevel.)) {
                  msLevel. <- base::sort(unique(msLevel(object)))
              } else {
                  if (!is.numeric(msLevel.))
                      stop("'msLevel' must be numeric!")
              }
              ps <- ProcessingStep("clean", list(all = all,
                                                 msLevel. = msLevel.))
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@cleaned <- TRUE
              object@processingData@processing <-
                  c(object@processingData@processing,
                    paste0("Spectra of MS level(s) ",
                           paste0(msLevel., sep = ", "),
                           " cleaned: ", date()))
              return(object)
          })

############################################################
## trimMz
##
## Add the "trimMz" ProcessingStep to the queue and update the
## processingData.
setMethod("trimMz", signature("OnDiskMSnExp", "numeric"),
          function(object, mzlim, msLevel.) {
              .Deprecated("filterMz")
              return(filterMz(object, mz = mzlim, msLevel.))
          })

############################################################
## normalize
##
## Handle the `normalize` method for OnMSnExp objects. We're adding a
## ProcessingStep for later, lazy processing. This will cause the
## `normalise` method to be applied to each spectrum once spectrum
## data (or intensity etc) is extracted.
setMethod("normalize", "OnDiskMSnExp",
          function(object, method = c("max", "sum"), ...) {
              method <- match.arg(method)
              ps <- ProcessingStep("normalise", list(method = method))
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@processing <-
                  c(object@processingData@processing,
                    paste0("Spectra normalised (", method, "): ",
                           date()))
              object@processingData@normalised <- TRUE
              return(object)
          })

############################################################
## bin
##
## That's a little tricky; we can't just add the `bin` Spectrum-method
## to the spectraProcessingQueue as this would be fairly inefficient
## (the breaks are calculated on the mz of the full data set).
## So:
## 1) need to get the mz range (from fData?); depending on the
##    msLevel.  argument however only from those spectra matching the
##    MS level.
## 2) add the mz range as a parameter to the processing queue.
## This is pretty slow, but should be robust. Eventually C-level
## binning might be faster.
setMethod("bin", "OnDiskMSnExp", function(x, binSize = 1L, msLevel.) {
    if (missing(msLevel.)) {
        msLevel. <- base::sort(unique(msLevel(x)))
    } else {
        if (!is.numeric(msLevel.))
            stop("'msLevel' must be numeric!")
    }
    ## Check if we have these MS levels
    if (!any(unique(msLevel(x)) %in% msLevel.)) {
        warning("No spectra of the specified MS level present.")
        return(x)
    }
    ## Get the M/Z range; note: calling spectrapply and returning just
    ## the M/Z range per spectrum is about twice as fast than getting
    ## all M/Z values and calculating the range on that (i.e.
    ## range(mz(object)))
    mzr <- range(unlist(spectrapply(filterMsLevel(x, msLevel. = msLevel.),
                                    FUN = function(z) {
                                        return(range(mz(z), na.rm = TRUE))
                                    }, BPPARAM = bpparam())))
    breaks <- seq(floor(mzr[1]), ceiling(mzr[2]), by = binSize)
    ## Now add the processing step
    ps <- ProcessingStep("bin", list(breaks = breaks, msLevel. = msLevel.))
    x@spectraProcessingQueue <- c(x@spectraProcessingQueue,
                                  list(ps))
    ## And add the processing info.
    x@processingData@processing <- c(x@processingData@processing,
                                     paste0("Spectra of MS level(s) ",
                                            paste0(msLevel., sep = ", "),
                                            " binned: ", date()))
    x
})

############################################################
## smooth
setMethod("smooth", "OnDiskMSnExp",
          function(x, method = c("SavitzkyGolay", "MovingAverage"),
                   halfWindowSize = 2L, verbose = isMSnbaseVerbose(),
                   msLevel. = unique(msLevel(x)), ...) {
              method <- match.arg(method)
              ps <- ProcessingStep("smooth",
                                   list(method = method,
                                        halfWindowSize = halfWindowSize,
                                        msLevel. = msLevel., ...))
              x@spectraProcessingQueue <- c(x@spectraProcessingQueue,
                                            list(ps))
              x@processingData@processing <- c(x@processingData@processing,
                                               paste0("Spectra smoothed (",
                                                      method, "): ",
                                                      date()))
              x@processingData@normalised <- TRUE
              return(x)
          })



############################################################
## pickPeaks
## Add a pickPeaks step to the processing queue
setMethod("pickPeaks", "OnDiskMSnExp",
          function(object, halfWindowSize = 3L,
                   method = c("MAD", "SuperSmoother"),
                   SNR = 0L, refineMz = c("none", "kNeighbors", "kNeighbours",
                                          "descendPeak"),
                   msLevel. = unique(msLevel(object)), ...) {
              method <- match.arg(method)
              refineMz <- match.arg(refineMz)
              ## If we're using an approach that combines spectra we have
              ## to add a specific processing step before.
              ps <- ProcessingStep("pickPeaks",
                                   list(method = method,
                                        halfWindowSize = halfWindowSize,
                                        SNR = SNR,
                                        ignoreCentroided = TRUE,
                                        refineMz = refineMz,
                                        msLevel. = msLevel., ...))
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object <- logging(
                  object, paste0("peak picking: ", method, " noise estimation",
                                 " and ", refineMz, " centroid m/z refinement",
                                 " on spectra of MS level(s)",
                                 paste0(msLevel., collapse = ", ")))
              object@processingData@smoothed <- TRUE
              fData(object)$centroided[msLevel(object) %in% msLevel.] <- TRUE
              if (validObject(object))
                  object
          })

############################################################
## compareSpectra
setMethod("compareSpectra", c("OnDiskMSnExp", "missing"),
          function(x, fun = c("common", "cor", "dotproduct"), ...) {
              fun <- match.arg(fun)
              ## res <- suppressMessages(compare_MSnExp(object1, fun, ...))
              ## return(res)
              ## Alternatively, we could get all of the spectra and doing
              ## the comparisons on the list; this might however turn out to
              ## be quite memory demanding...
              sps <- spectra(x)
              nm <- featureNames(x)
              cb <- combn(nm, 2, function(z) {
                  compare_Spectra(sps[[z[1]]], sps[[z[2]]], fun = fun, ...)
              })
              m <- matrix(NA, length(x), length(x),
                          dimnames = list(nm, nm))
              ## fill lower triangle of the matrix
              m[lower.tri(m)] <- cb
              ## copy to upper triangle
              for (i in 1:nrow(m)) {
                  m[i, ] <- m[, i]
              }

              return(m)
          })

############################################################
## estimateNoise
##
## estimateNoise directly estimates the noise and returns it. We're going to
## call the method using the spectrapply method.
## Note: would also work directly with the MSnExp implementation, but this would
## require that all spectra are loaded first, while here we ensure parallel
## processing by file; should be faster and memory efficient.
setMethod("estimateNoise", "OnDiskMSnExp",
          function(object, method = c("MAD", "SuperSmoother"), ...) {
              method <- match.arg(method)
              ## Call estimateNoise_Spectrum using the spectrapply.
              res <- spectrapply(object, FUN = estimateNoise_Spectrum,
                                 method = method, ...)
              return(res)
          })

############################################################
## removeReporters
##
## Adds the removeReporters Spectrum processing step to the lazy processing
## queue
setMethod("removeReporters", "OnDiskMSnExp",
          function(object, reporters = NULL, clean = FALSE,
                   verbose = isMSnbaseVerbose()) {
              if (is.null(reporters)) {
                  return(object)
              } else {
                  ## Check if we've got msLevel > 1.
                  if (all(msLevel(object) == 1))
                      stop("No MS level > 1 spectra present! Reporters can",
                           " only be removed from spectra >= 2!")
                  ps <- ProcessingStep("removeReporters",
                                       list(reporters = reporters,
                                            clean = clean,
                                            suppressWarnings = TRUE))
                  object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                     list(ps))
                  repname <- names(reporters)
                  object@processingData@processing <- c(object@processingData@processing,
                                                   paste("Removed",
                                                         repname,
                                                         "reporter ions ",
                                                         date()))
                  return(object)
              }
          })

############################################################
## spectrapply
##
## That's the main method to apply functions to the object's spectra, or
## to just return a list with the spectra, if FUN is empty.
## Parallel processing by file can be enabled using BPPARAM.
setMethod("spectrapply", "OnDiskMSnExp", function(object, FUN = NULL,
                                                  BPPARAM = bpparam(), ...) {
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)
    ## Get the fastLoad option.
    fast_load <- isMSnbaseFastLoad()
    isOK <- validateFeatureDataForOnDiskMSnExp(fData(object))
    if (!is.null(isOK))
        stop(isOK)
    fDataPerFile <- split.data.frame(fData(object),
                                     f = fData(object)$fileIdx)
    fNames <- fileNames(object)
    theQ <- processingQueue(object)
    vals <- bplapply(fDataPerFile,
                     FUN = .applyFun2SpectraOfFileMulti,
                     filenames = fNames,
                     queue = theQ,
                     APPLYFUN = FUN,
                     fastLoad = fast_load,
                     BPPARAM = BPPARAM,
                     ...)
    names(vals) <- NULL
    vals <- unlist(vals, recursive = FALSE)
    vals[rownames(fData(object))]
})

setMethod("isolationWindowLowerMz", "OnDiskMSnExp", function(object) {
    if (all(c("isolationWindowTargetMZ", "isolationWindowLowerOffset") %in%
            colnames(fData(object))))
        return(fData(object)$isolationWindowTargetMZ -
                            fData(object)$isolationWindowLowerOffset)
    rep(NA_real_, length(object))
})
setMethod("isolationWindowUpperMz", "OnDiskMSnExp", function(object) {
    if (all(c("isolationWindowTargetMZ", "isolationWindowUpperOffset") %in%
            colnames(fData(object))))
        return(fData(object)$isolationWindowTargetMZ +
                            fData(object)$isolationWindowUpperOffset)
    rep(NA_real_, length(object))
})

setMethod("splitByFile", c("OnDiskMSnExp", "factor"), function(object, f) {
    .on_disk_split_by_file(object)
})
