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
          lockEnvironment(ans@assayData, bindings = TRUE)
          ans <- MSnbase:::logging(ans, "Coerced from OnDiskMSnExp")
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
          function(object) {
              val <- fData(object)$centroided
              names(val) <- featureNames(object)
              return(val)
          })

setReplaceMethod("centroided",
                 signature(object = "OnDiskMSnExp", value = "logical"),
                 function(object, value, msLevel.) {
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

## setMethod("isCentroided", "OnDiskMSnExp",
##           function(object, ..., verbose = TRUE) {
##               pkl <- lapply(spectra(object), as.data.frame)
##               ctrd <- lapply(pkl, .isCentroided, ...)
##               ctrd <- unlist(ctrd, use.names = FALSE)
##               if (verbose) print(table(ctrd, msLevel(object)))
##               ctrd
##           })

setMethod("isCentroided", "OnDiskMSnExp",
          function(object, ..., verbose = TRUE) {
              res <- spectrapply(object, function(z, ...) {
                  return(.isCentroided(as(z, "data.frame"), ...))
              }, ...)
              ctrd <- unlist(res, use.names = FALSE)
              if (verbose) print(table(ctrd, msLevel(object)))
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
                 function(object, value, msLevel.) {
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
          function(object, BPPARAM = bpparam()) {
    skipFun <- c("pickPeaks")
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
    if (recalc) {
        vals <- unlist(spectrapply(object,
                                   FUN = tic,
                                   BPPARAM = BPPARAM))
    } else {
        vals <- fData(object)$totIonCurrent
    }
    names(vals) <- featureNames(object)
    return(vals)
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
              if (is.character(i))
                  i <- match(i, featureNames(x))
              if (any(is.na(i)))
                  stop("subscript out of bounds")
              ## Keep the fromFile information; that would be
              ## changed by [ subsetting.
              fromF <- fromFile(x)[i]
              ## Subsetting first the OnDiskMSnExp with subsequent
              ## spectra extraction is considerably faster.
              spctr <- spectra(x[i])
              ## Re-setting the original fromFile values.
              spctr <- mapply(spctr, fromF, FUN = function(y, z) {
                  y@fromFile <- z
                  return(y)
              })
              if (length(spctr) == 1)
                  spctr <- spctr[[1]]
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
                   verbose = TRUE,
                   ...) {
              method <- match.arg(method)
              if (!all(msLevel(object) == 2)) {
                  message("Currently only MS2 quantitation: filtering MS2 spectra.")
                  object <- filterMsLevel(object, msLevel. = 2L)
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
setMethod("[", signature(x = "OnDiskMSnExp",
                         i = "logicalOrNumeric",
                         j = "missing",
                         drop = "missing"),
          function(x, i, j, drop) {
              if (!(is.logical(i) | is.numeric(i)))
                  stop("Subsetting works only with numeric or logical!")
              if (is.logical(i)) {
                  if (length(i) != nrow(fData(x)))
                      stop("If 'i' is logical its length has to match",
                           " the number of spectra!")
                  i <- which(i)
              }
              i <- sort(i)  ## Force sorted!
              ## Now subset the featureData. The function will
              ## complain if i is outside the range.
              if (length(i) == 0) x@featureData <- x@featureData[0, ]
              else x@featureData <- subsetFeatureDataBy(featureData(x),
                                                        index = i)
              ## Check also that processingData and experimentData
              ## match the *new* featureData:
              file <- sort(unique(fromFile(x)))
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
              expD@instrumentManufacturer <- expD@instrumentManufacturer[file]
              expD@instrumentModel <- expD@instrumentModel[file]
              expD@ionSource <- expD@ionSource[file]
              expD@analyser <- expD@analyser[file]
              expD@detectorType <- expD@detectorType[file]
              x@experimentData <- expD
              ## Update fromFile in spectra/featureData.
              fromFile(x) <- match(fromFile(x), file)
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
              object <- object[sort(gotEm)]
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
          function(object, t = "min", verbose = TRUE, msLevel.) {
              if (missing(t))
                  t <- "min"
              if (!is.numeric(t)) {
                  if (t != "min")
                      stop("Argument 't' has to be either numeric or 'min'!")
              }
              if (missing(msLevel.)) {
                  msLevel. <- sort(unique(msLevel(object)))
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
          function(object, all = FALSE, verbose = TRUE, msLevel.) {
              if (!is.logical(all))
                  stop("Argument 'all' is supposed to be a logical!")
              if (missing(msLevel.)) {
                  msLevel. <- sort(unique(msLevel(object)))
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
setMethod("bin", "OnDiskMSnExp", function(object, binSize = 1L, msLevel.) {
    if (missing(msLevel.)) {
        msLevel. <- sort(unique(msLevel(object)))
    } else {
        if (!is.numeric(msLevel.))
            stop("'msLevel' must be numeric!")
    }
    ## Check if we have these MS levels
    if (!any(unique(msLevel(object)) %in% msLevel.)) {
        warning("No spectra of the specified MS level present.")
        return(object)
    }
    ## Get the M/Z range; note: calling spectrapply and returning just
    ## the M/Z range per spectrum is about twice as fast than getting
    ## all M/Z values and calculating the range on that (i.e.
    ## range(mz(object)))
    mzr <- range(unlist(spectrapply(filterMsLevel(object, msLevel. = msLevel.),
                                    FUN = function(z) {
                                        return(range(mz(z), na.rm = TRUE))
                                    }, BPPARAM = bpparam())))
    breaks <- seq(floor(mzr[1]), ceiling(mzr[2]), by = binSize)
    ## Now add the processing step
    ps <- ProcessingStep("bin", list(breaks = breaks, msLevel. = msLevel.))
    object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                       list(ps))
    ## And add the processing info.
    object@processingData@processing <- c(object@processingData@processing,
                                          paste0("Spectra of MS level(s) ",
                                                 paste0(msLevel., sep = ", "),
                                                 " binned: ", date()))
    return(object)
})

############################################################
## smooth
setMethod("smooth", "OnDiskMSnExp",
          function(x, method = c("SavitzkyGolay", "MovingAverage"),
                   halfWindowSize = 2L, verbose = TRUE, ...) {
              method <- match.arg(method)
              ps <- ProcessingStep("smooth",
                                   list(method = method,
                                        halfWindowSize = halfWindowSize, ...))
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
                   SNR = 0L, ...) {
              method <- match.arg(method)
              ps <- ProcessingStep("pickPeaks",
                                   list(method = method,
                                        halfWindowSize = halfWindowSize,
                                        SNR = SNR,
                                        ignoreCentroided = TRUE, ...))
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              object@processingData@processing <-
                  c(object@processingData@processing,
                    paste0("Peak picking (", method, "): ", date()))
              object@processingData@smoothed <- TRUE
              fData(object)$centroided <- TRUE
              if (validObject(object))
                  return(object)
          })

############################################################
## compareSpectra
setMethod("compareSpectra", c("OnDiskMSnExp", "missing"),
          function(object1, fun = c("common", "cor", "dotproduct"), ...) {
              fun <- match.arg(fun)
              ## res <- suppressMessages(compare_MSnExp(object1, fun, ...))
              ## return(res)
              ## Alternatively, we could get all of the spectra and doing
              ## the comparisons on the list; this might however turn out to
              ## be quite memory demanding...
              sps <- spectra(object1)
              nm <- featureNames(object1)
              cb <- combn(nm, 2, function(x) {
                  compare_Spectra(sps[[x[1]]], sps[[x[2]]], fun = fun, ...)
              })
              m <- matrix(NA, length(object1), length(object1),
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
## require that all spectra are loaded first, while here we ensure parallel processin
## by file; should be faster and memory efficient.
setMethod("estimateNoise", "OnDiskMSnExp",
          function(object, method = c("MAD", "SuperSmoother"), ...) {
              method <- match.arg(method)
              ## Call estimateNoise_Spectrum using the spectrapply.
              res <- spectrapply(object, FUN = estimateNoise_Spectrum,
                                 method = method, ...)
              return(res)
          })



##============================================================
##  --  HELPER FUNCTIONS  --
##
##------------------------------------------------------------


## Apply a function to the spectra of a file.
## The function will first read the raw data, create Spectrum objects
## from it, apply all ProcessingSteps apply the specified function and
## return its results. This is done on a per-file basis, so we can do
## that in parallel.
## Input args:
## fData: the data.frame containing all feature data info for the
##   current file. All spectra specified by column "spIdx" are read.
## filenames: the files from which the data should be read; this
##   should be all file names (i.e. those returned) by
##   fileNames(object), as we're using the fileIdx column in fData to
##   get the "real" filename.
## queue: a list of ProcessingStep objects that should be applied to
##   the Spectrum.
## APPLYFUN: the function that should be applied; if APPLYFUN is NULL
##   it returns the spectra!
.applyFun2SpectraOfFile <- function(fData, filenames,
                                    queue = NULL,
                                    APPLYFUN = NULL) {
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are requierd!")
    ## Subset the fData based on scnIDx; if provided.
    filename <- filenames[fData[1, "fileIdx"]]
    if (any(fData$msLevel > 1))
        stop("on-the-fly import currently only supported for MS1 level data.")
    ## Open the file.
    fileh <- mzR::openMSfile(filename)
    on.exit(expr = mzR::close(fileh))
    ## Now run the stuff per spectrum, i.e. read data, create object,
    ## apply the fun.  Note that we're splitting the matrix not the
    ## data.frame, as that's faster.
    res <-
        lapply(split(fData[, c("fileIdx", "spIdx", "centroided",
                               "acquisitionNum", "peaksCount",
                               "totIonCurrent", "retentionTime",
                               "polarity")],
                     rownames(fData)),
               FUN = function(z, theQ, fh, APPLYFUN){
                   ## Read the data.
                   spD <- mzR::peaks(fh, z[1, 2])
                   sp <- Spectrum1(peaksCount = z[1, 5],
                                   rt = z[1, 7],
                                   acquisitionNum = z[1, 4],
                                   scanIndex = z[1, 2],
                                   mz = spD[, 1],
                                   intensity = spD[, 2],
                                   fromFile = z[1, 1],
                                   centroided = z[1, 3],
                                   polarity = z[1, 8])
                   ## Now, apply the Queue.
                   if (length(theQ) > 0){
                       for (pStep in theQ){
                           sp <- execute(pStep, sp)
                       }
                   }
                   if (is.null(APPLYFUN))
                       return(sp)
                   ## Now what remains is to apply the APPLYFUN and
                   ## return results.
                   return(APPLYFUN(sp))
               }, theQ = queue, fh = fileh, APPLYFUN = APPLYFUN)
    return(res)
}
## The same function but using the standard "new" constructor from R.
.applyFun2SpectraOfFileSlow <- function(fData, filenames,
                                        queue = NULL, APPLYFUN = NULL){
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are requierd!")
    filename <- filenames[fData[1, "fileIdx"]]
    ## Open the file.
    fileh <- mzR::openMSfile(filename)
    on.exit(expr = mzR::close(fileh))

    msLevel1 <- which(fData$msLevel == 1)
    msLevelN <- which(fData$msLevel > 1)
    ## Process MS1 and MSn separately
    if (length(msLevel1) >= 1) {
        ms1fd <- fData[msLevel1, , drop = FALSE]
        ## Now run the stuff per spectrum, i.e. read data, create
        ## object, apply the fun.  Note that we're splitting the
        ## matrix not the data.frame, as that's faster.
        res <- lapply(split(ms1fd, rownames(ms1fd)),
                      FUN = function(z, theQ, fh, APPLYFUN) {
                          ## Read the data.  According to issue #103
                          ## we should use acquisitionNum, not
                          ## spectrum idx.
                          spD <- mzR::peaks(fh, z$acquisitionNum)
                          sp <- new("Spectrum1",
                                    fromFile = z$fileIdx,
                                    scanIndex = z$spIdx,
                                    centroided = z$centroided,
                                    acquisitionNum = z$acquisitionNum,
                                    ## peaksCount = z[1,5],
                                    rt = z$retentionTime,
                                    polarity = z$polarity,
                                    mz = spD[, 1],
                                    intensity = spD[, 2])
                          ## Now, apply the Queue.
                          if (length(theQ) > 0) {
                              for (i in 1:length(theQ)) {
                                  sp <- execute(theQ[[i]], sp)
                              }
                          }
                          if (is.null(APPLYFUN))
                              return(sp)
                          ## Now what remains is to apply the APPLYFUN
                          ## and return results.
                          return(APPLYFUN(sp))
                      }, theQ = queue, fh = fileh, APPLYFUN = APPLYFUN)
    } else {
        res <- list()
    }
    if (length(msLevelN) >= 1) {
        msnfd <- fData[msLevelN, , drop = FALSE]
        ## TODO: write/use C-constructor
        ## For now we're using the lapply, new() approach iteratively
        ## reading each spectrum from file and creating the Spectrum2.
        res2 <- lapply(split(msnfd, rownames(msnfd)),
                       function(z, theQ, fh, APPLYFUN) {
                           spectD <- mzR::peaks(fh, z$acquisitionNum)
                           sp <- new("Spectrum2",
                                     merged = z$mergedScan,
                                     precScanNum = z$precursorScanNum,
                                     precursorMz = z$precursorMZ,
                                     precursorIntensity = z$precursorIntensity,
                                     precursorCharge = z$precursorCharge,
                                     collisionEnergy = z$collisionEnergy,
                                     msLevel = z$msLevel,
                                     rt = z$retentionTime,
                                     acquisitionNum = z$acquisitionNum,
                                     scanIndex = z$spIdx,
                                     mz = spectD[, 1],
                                     intensity = spectD[, 2],
                                     fromFile = z$fileIdx,
                                     centroided = z$centroided,
                                     polarity = z$polarity)
                           ## Now, apply the Queue.
                           if (length(theQ) > 0) {
                               for (i in 1:length(theQ)) {
                                   sp <- execute(theQ[[i]], sp)
                               }
                           }
                           if (is.null(APPLYFUN))
                               return(sp)
                           ## Now what remains is to apply the
                           ## APPLYFUN and return results.
                           return(APPLYFUN(sp))
                       }, theQ = queue, fh = fileh, APPLYFUN = APPLYFUN)
        names(res2) <- rownames(fData)[msLevelN]
        res <- c(res, res2)
    }
    ## Ensure that ordering is the same than in fData:
    res <- res[match(rownames(fData), names(res))]
    return(res)
}
## Using the C constructor that takes all values at once and creates a
## list of Spectrum1 objects, applies processing steps, applies the
## provided function and returns its results - or the list of
## Spectrum1 objects if APPLYFUN = NULL.
## Arguments:
## o fData: either a full data.frame (returned by fData(OnDiskMSnExp))
##   or a sub-set forspecific spectra. The data.frame should ONLY
##   CONTAIN VALUES FOR SPECTRA OF ONE FILE!
## o filenames: fileNames(object)
## o queue: object@spectraProcessingQueue; if lenght > 0 all
##   processing steps will be applied to the created Spectrum1
##   objects.
## o APPLYFUN: the function to be applied to the Spectrum1 objects
##   (such as ionCount etc).  If NULL the function returns the list of
##   Spectrum1 objects.
## o APPLYARGS: additional arguments for the APPLYFUN
.applyFun2SpectraOfFileMulti <- function(fData, filenames,
                                         queue = NULL,
                                         APPLYFUN = NULL,
                                         ...) {
    if (missing(fData) | missing(filenames))
        stop("Both 'fData' and 'filenames' are required!")
    filename <- filenames[fData[1, "fileIdx"]]
    ## if (any(fData$msLevel > 1))
    ##     stop("on-the-fly import currently only supported for MS1 level data.")
    ## Open the file.
    fileh <- mzR::openMSfile(filename)
    hd <- header(fileh)
    on.exit(expr = mzR::close(fileh))
    msLevel1 <- which(fData$msLevel == 1)
    msLevelN <- which(fData$msLevel > 1)
    ## Process MS1 and MSn separately
    if (length(msLevel1) >= 1) {
        ms1fd <- fData[msLevel1, , drop = FALSE]
        ## Reading all of the data in "one go". According to issue
        ## #103 we should use acquisitionNum, not spectrum idx.
        ## See issue #118 for an explanation of the match
        allSpect <- mzR::peaks(fileh,
                               match(ms1fd$acquisitionNum, hd$acquisitionNum))
        ## If we have more than one spectrum the peaks function returns a list.
        if (is(allSpect, "list")) {
            nValues <- lengths(allSpect) / 2
            allSpect <- do.call(rbind, allSpect)
        } else {
            ## otherwise it's a matrix, e.g. if only a single scan
            ## index was provided.
            nValues <- nrow(allSpect)
        }
        ## Call the C-constructor to create a list of Spectrum1
        ## objects.
        res <- Spectra1(peaksCount = nValues,
                        scanIndex = ms1fd$spIdx,
                        rt = ms1fd$retentionTime,
                        acquisitionNum = ms1fd$acquisitionNum,
                        mz = allSpect[, 1],
                        intensity = allSpect[, 2],
                        centroided = ms1fd$centroided,
                        smoothed = ms1fd$smoothed,
                        fromFile = ms1fd$fileIdx,
                        polarity = ms1fd$polarity,
                        tic = ms1fd$totIonCurrent,
                        nvalues = nValues)
        names(res) <- rownames(ms1fd)
    } else {
        res <- list()
    }
    if (length(msLevelN) >= 1) {
        msnfd <- fData[msLevelN, , drop = FALSE]
        ## Reading all of the data in "one go".
        ## See issue #118 for an explanation of the match
        allSpect <- mzR::peaks(fileh,
                               match(msnfd$acquisitionNum, hd$acquisitionNum))
        ## If we have more than one spectrum the peaks function returns a list.
        if (is(allSpect, "list")) {
            nValues <- lengths(allSpect) / 2
            allSpect <- do.call(rbind, allSpect)
        } else {
            ## otherwise it's a matrix, e.g. if only a single scan
            ## index was provided.
            nValues <- nrow(allSpect)
        }
        ## Call the C-constructor to create a list of Spectrum2
        ## objects.
        res2 <- Spectra2(peaksCount = nValues,
                         scanIndex = msnfd$spIdx,
                         tic = msnfd$totIonCurrent,
                         rt = msnfd$retentionTime,
                         acquisitionNum = msnfd$acquisitionNum,
                         mz = allSpect[, 1],
                         intensity = allSpect[, 2],
                         centroided = msnfd$centroided,
                         fromFile = msnfd$fileIdx,
                         polarity = msnfd$polarity,
                         nvalues = nValues,
                         msLevel = msnfd$msLevel,
                         merged = msnfd$mergedScan,
                         precScanNum = msnfd$precursorScanNum,
                         precursorMz = msnfd$precursorMZ,
                         precursorIntensity = msnfd$precursorIntensity,
                         precursorCharge = msnfd$precursorCharge,
                         collisionEnergy = msnfd$collisionEnergy)
        names(res2) <- rownames(msnfd)
        res <- c(res, res2)
    }
    ## Ensure that ordering is the same than in fData:
    res <- res[match(rownames(fData), names(res))]
    ## If we have a non-empty queue, we might want to execute that too.
    if (!is.null(APPLYFUN) | length(queue) > 0){
        if (length(queue) > 0) {
            message("Apply lazy processing step(s):")
            for (j in 1:length(queue))
                message(" o '", queue[[j]]@FUN, "' with ",
                        length(queue[[j]]@ARGS), " argument(s).")
        }
        res <- lapply(res, function(z, theQ, APPLF, ...){
            if (length(theQ) > 0) {
                for (pStep in theQ) {
                    z <- execute(pStep, z)
                }
            }
            if (is.null(APPLF))
                return(z)
            ## return(APPLF(z, APPLA))
            return(do.call(APPLF, args = c(list(z), ...)))
        }, theQ = queue, APPLF = APPLYFUN, ...)
    }
    return(res)
}
