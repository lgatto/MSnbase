## setMethod("initialize",
##           signature(.Object="pSet"),
##           function(.Object, ..., .cache) {
##               cat("initialize pSet")
##               if (missing(.cache)) {
##                   .cache <- new.env()
##                   assign("level", 0 ,.cache)
##                   lockEnvironment(.cache)
##               }
##               callNextMethod(.Object, ..., .cache = .cache)
##           })

## That (according to the R-help) should be the "correct" way to set specific
## slot values if not provided. With the "old" version above some slot values (like
## 'onDisk') were not propagated.
setMethod("initialize",
          "pSet",
          function(.Object, ...) {
              if (!any(names(list(...)) == ".cache")) {
                  .cache <- new.env()
                  assign("level", 0, .cache)
                  lockEnvironment(.cache)
                  .Object@.cache <- .cache
              }
              callNextMethod()
          })

setValidity("pSet", function(object) {
    msg <- validMsg(NULL, NULL)
    ## Skip some (most) of the tests for a OnDiskMSnExp, since we don't have
    ## spectrum data available, i.e. assayData is empty.
    if (length(object)) {
        if (!all(sapply(assayData(object), function(x) inherits(x, "Spectrum"))))
            msg <- validMsg(msg,
                            "assayData must contain 'Spectrum' objects.")
        msl <- msLevel(object)
        if (length(unique(msl)) > 1)
            warning(paste0("Different MS levels in ", class(object),
                           " object: ", unique(msl)))
        ## checking number of spectra in assayData and
        ##          number of features in featureData
        nspectra  <- length(assayData(object))
        nfeatures <- nrow(featureData(object))
        if (nspectra != nfeatures)
            msg <- validMsg(msg, "Unequal number of spectra in assayData and features in featureData.")
        if (length(spectra(object)) != length(ls(assayData(object))))
            msg <- validMsg(msg, "Object size inconsistence using assayData() and spectra() methods.")
        if (!identical(featureNames(object),
                       ls(assayData(object))))
            msg <- validMsg(msg, "featureNames differ between assayData and featureData.")
        ## checking number of files in phenoData and
        ##          number of files in assayData
        ## removing the following check because we add identification files
        ## (a modified version using "<" instead of "!=" is below
        ## nfilesprocData   <- length(processingData(object)@files)
        ## nfilesSpectra <- length(unique(unlist(eapply(assayData(object),fromFile))))
        ## if (nfilesprocData != nfilesSpectra)
        ##   msg <- validMsg(msg, "unequal number of files in assayData and processingData.")
        aFileIds <- fromFile(object)
        fFileIds <- fData(object)$file
        if (length(fFileIds) && any(aFileIds != fFileIds))
            msg <- validMsg(msg, "Mismatch of files in assayData and processingData.")
        nfilesprocData   <- length(processingData(object)@files)
        nfilesSpectra <- length(unique(aFileIds))
        if (nfilesprocData < nfilesSpectra)
            msg <- validMsg(msg, "More spectra files in assayData than in processingData.")
        if (length(sampleNames(object)) != nrow(pData(object)))
            msg <- validMsg(msg, "Different number of samples accoring to sampleNames and pData.")
        ## phenoData: at this stage, we don't know how many sample there will be
        ## expecting one for label-free, but we could have 4, 6, 8 for isobaric tagging
        ## or 2 for metabolic labelling
        ## if (nrow(pData(object)) > 1)
        ##     message("Detected ", nrow(pData(object)), " samples in pData().")
        ##
        ## protocolData not checked yet - depends very much
        ## on type of assay (MS1, MS2 quant, reporter ions, ...)
        if (!cacheEnvIsLocked(object))
            msg <- validMsg(msg, "'.cache' environment is not locked.")
        if (!exists("level", envir = object@.cache))
            msg <- validMsg(msg, "'.cache' level not defined.")
        hd <- header(object)
        if (nrow(hd) != length(object))
            msg("(Cached) header nrow and object length differ.")
        sapply(spectra(object), validObject)
    }
    if (is.null(msg)) TRUE else msg
})


setMethod("[", "pSet",
          function(x, i, j = "missing", drop = "missing") {
              if (!(is.logical(i) | is.numeric(i)))
                  stop("subsetting works only with numeric or logical")
              if (is.numeric(i)) {
                  if (max(i) > length(x) | min(i) < 1)
                      stop("subscript out of bounds")
                  i <- base::sort(i) ## crash if unsorted (because of
                  ## (alphanumerical) order in ls(assayData(.))  and
                  ## featureNames(featureData) that have to be
                  ## indentical - see issues #70 and #71)
              }
              whichElements <- ls(assayData(x))[i]
              sel <- featureNames(x) %in% whichElements
              ## Sub-setting also processingData, phenoData and experimentData.
              file <- base::sort(unique(fromFile(x)[sel]))
              ## o Sub-set the phenoData.
              pd <- phenoData(x)[file, , drop = FALSE]
              pData(pd) <- droplevels(pData(pd))
              x@phenoData <- pd
              ## Sub-set the files.
              ## This one should never be empty
              x@processingData@files <- x@processingData@files[file]
              ## Sub-set the experimentData: these slots can be empty
              ## character(), producing NA if subest with character()[file]
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
              ## Fix the fromFile property
              newFromFile <- base::match(fromFile(x), file)
              names(newFromFile) <- names(fromFile(x))
              ## Proceed.
              orghd <- header(x)
              olde <- assayData(x)
              newe <- new.env(parent = emptyenv())
              if (length(whichElements) > 0) {
                  for (el in whichElements) {
                      sp <- olde[[el]]
                      sp@fromFile <- unname(newFromFile[el])
                      newe[[el]] <- sp
                  }
                  if (environmentIsLocked(olde))
                      lockEnvironment(newe,
                                      bindings = bindingIsLocked(el, olde))
              } else {
                      lockEnvironment(newe, bindings = TRUE)
              }
              x@assayData <- newe
              x@featureData <- featureData(x)[i, ]
              if (is.logical(i)) {
                  x@processingData@processing <-
                      c(processingData(x)@processing,
                        paste("Data [logically] subsetted ", sum(i),
                              " spectra: ", date(), sep = ""))
              } else if (is.numeric(i)) {
                  x@processingData@processing <-
                      c(processingData(x)@processing,
                        paste("Data [numerically] subsetted ",
                              length(i), " spectra: ", date(),
                              sep = ""))
              } else {
                  x@processingData@processing <-
                      c(processingData(x)@processing,
                        paste("Data subsetted ", i ,": ", date(),
                              sep = ""))
              }
              if (x@.cache$level > 0) {
                  ## no caching for empty objects
                  .cache <- ifelse(length(x) > 1, x@.cache$level, 0)
                  x@.cache <- setCacheEnv(list(assaydata = assayData(x),
                                               hd = orghd[sel, ]), ## faster than .header for big instances
                                          .cache)
              }
              if (validObject(x))
                  return(x)
          })


setMethod("[[", "pSet",
          function(x, i, j = "missing", drop = "missing") {
              if (length(i) != 1)
                  stop("subscript out of bounds")
              if (!is.character(i))
                  i <- featureNames(x)[i]
              return(get(i, envir = assayData(x)))
          })

setMethod("precursorMz", "pSet",
          function(object) {
              ## this assumes that if first spectrum
              ## has msLevel > 1, all have
              if (.firstMsLevel(object) > 1)
                  return(sapply(spectra(object), precursorMz))
              stop("No precursor MZ value for MS1 spectra.")
          })

## a general version
## setMethod("precursorMz","pSet",
##           function(object) {
##               msl <- msLevel(object)
##               ans <- rep(NA_integer_, length(object))
##               obj2 <- object[msl > 1]
##               ans2 <- sapply(spectra(obj2), slot, "precursorMz")
##               ans[msl > 1] <- ans2
##               names(ans) <- featureNames(object)
##               ans
##           })

setMethod("precScanNum", "pSet",
          function(object) {
              if (.firstMsLevel(object) > 1)
                  return(unlist(sapply(spectra(object), precScanNum)))
              stop("This experiment contains MS1 spectra.")
          })


setMethod("tic", "pSet",
          function(object) sapply(spectra(object), tic))

setMethod("ionCount", "pSet",
          function(object) sapply(spectra(object), ionCount))

setMethod("precursorCharge", "pSet",
          function(object) {
              if (.firstMsLevel(object) > 1)
                  return(sapply(spectra(object), precursorCharge))
              stop("No precursor MZ value for MS1 spectra.")
          })

setMethod("precursorIntensity", "pSet",
          function(object) {
              if (.firstMsLevel(object) > 1)
                  return(sapply(spectra(object), precursorIntensity))
              stop("No precursor data for MS1 spectra.")
          })

setMethod("acquisitionNum", "pSet",
          function(object) sapply(spectra(object), acquisitionNum))
setMethod("scanIndex", "pSet",
          function(object) sapply(spectra(object), scanIndex))

setMethod("rtime", "pSet",
          function(object) sapply(spectra(object), rtime))

setMethod("centroided", "pSet",
          function(object, na.fail = FALSE)
              sapply(spectra(object), centroided, na.fail))

setReplaceMethod("centroided",
                 signature(object = "pSet",
                           value = "logical"),
                 function(object, value) {
                     if (length(value) == 1)
                         value <- rep(value, length(object))
                     if (length(object) != length(value))
                         stop("Length of replacement value is different than number of spectra.")
                     sl <- spectra(object)
                     for (i in 1:length(sl))
                         centroided(sl[[i]]) <- value[i]
                     object@assayData <- as.environment(sl)
                     if (validObject(object))
                         return(object)
                 })

setMethod("smoothed", "pSet",
          function(object) sapply(spectra(object), smoothed))

setReplaceMethod("smoothed",
                 signature(object = "pSet",
                           value = "logical"),
                 function(object, value) {
                     if (length(value) == 1)
                         value <- rep(value, length(object))
                     if (length(object) != length(value))
                         stop("Length of replacement value is different than number of spectra.")
                     sl <- spectra(object)
                     for (i in 1:length(sl))
                         smoothed(sl[[i]]) <- value[i]
                     object@assayData <- as.environment(sl)
                     if (validObject(object))
                         return(object)
                 })

setMethod("peaksCount",
          signature("pSet", "missing"),
          function(object) sapply(spectra(object), peaksCount))

setMethod("peaksCount",
          signature("pSet", "numeric"),
          function(object, scans) {
              if (length(scans) == 1)
                  return(peaksCount(object[[scans]]))
              sapply(spectra(object)[scans], peaksCount)
          })

setMethod("msLevel", "pSet",
          function(object) sapply(spectra(object), msLevel))

setMethod("collisionEnergy", "pSet",
          function(object) {
              if (.firstMsLevel(object) > 1)
                  return(sapply(spectra(object), collisionEnergy))
              stop("No collision energy for MS1 spectra.")
          })

setMethod("intensity", "pSet",
          function(object) lapply(spectra(object),intensity))

setMethod("mz", "pSet",
          function(object) lapply(spectra(object),mz))

setMethod("polarity", "pSet",
          function(object) sapply(spectra(object), polarity))


setMethod("fromFile", "pSet",
          function(object) return(sapply(spectra(object), fromFile)))
setReplaceMethod("fromFile",
                 signature(object = "pSet",
                           value = "integer"),
                 function(object, value) {
                     if (length(object) != length(value))
                         stop("Length of replacement value is different than number of spectra.")
                     sl <- spectra(object)
                     for (i in 1:length(sl))
                         sl[[i]]@fromFile <- value[i]
                     object@assayData <- as.environment(sl)
                     if (validObject(object))
                         return(object)
                 })

setMethod("header",
          signature("pSet", "missing"),
          function(object) {
              ifelse(object@.cache$level > 0,
                     hd <- object@.cache$hd,
                     hd <- .header(object))
              return(hd)
          })

setMethod("header",
          signature = c("pSet", "numeric"),
          function(object, scans) {
              hd <- header(object)
              return(hd[scans, ])
          })

##################################################################
## NOTE: pSet is based in the eSet class defined in the Biobase
##       package, and accessor and replace methods for the common
##       slots are expected to have identical behaviours.
##       In these cases, code is heavily inspired and sometimes
##       directly copies from the respective eSet method.


## returns the dimensions of PhenoData
setMethod("dim", "pSet", function(x) dim(pData(x)))

## returns the number of spectra in the AssayData env
setMethod("length", "pSet", function(x) length(assayData(x)))

setMethod("assayData", "pSet", function(object) object@assayData)

setMethod("spectra", "MSnExp", function(object, ...) {
    sl <- as.list(assayData(object))
    fnames <- featureNames(object)
    ## reordering the spectra in the spectra list to match
    ## their order in featureData
    return(sl[fnames])
})

## setReplaceMethod("assayData",
##                  signature=signature(
##                    object="pSet",
##                    value="AssayData"),
##                  function(object, value) {
##                    object@assayData <- value ## may be lock env?
##                    return(object)
##                  })


setMethod("sampleNames",
          signature(object = "pSet"),
          function(object) sampleNames(phenoData(object)))

setMethod("fileNames",
          signature(object = "pSet"),
          function(object) processingData(object)@files)

## setReplaceMethod("fileNames",
##           signature(object="pSet", value="character"),
##           function(object, value) {
##             fileNames(object@processingData) <- value
##             return(object)
##           })

setReplaceMethod("sampleNames",
                 signature = signature(object = "pSet", value = "character"),
                 function(object, value) {
                     pd <- phenoData(object)
                     sampleNames(pd) <- value
                     prd <- protocolData(object)
                     if (nrow(prd) == 0) {
                         prd <- pd[, integer(0)]
                     } else {
                         sampleNames(prd) <- value
                     }
                     object@phenoData <- pd
                     object@protocolData <- prd
                     if (validObject(object))
                         return(object)
                 })

setMethod("featureNames",
          signature = signature(object = "pSet"),
          function(object) featureNames(featureData(object)))

setMethod("phenoData", "pSet", function(object) object@phenoData)
setMethod("pData", "pSet", function(object) pData(phenoData(object)))
setMethod("varMetadata",
          signature = signature(object = "pSet"),
          function(object) varMetadata(phenoData(object)))
setMethod("varLabels",
          signature = signature(object = "pSet"),
          function(object) varLabels(phenoData(object)))
setMethod("featureData",
          signature(object = "pSet"),
          function(object) object@featureData)

setReplaceMethod("featureData",
                 signature = signature(
                     object = "pSet",
                     value = "AnnotatedDataFrame"),
                 function(object, value) {
                     object@featureData <- value
                     if (validObject(object))
                         return(object)
                 })

setMethod("fData",
          signature = signature(object = "pSet"),
          function(object) pData(featureData(object)))

setReplaceMethod("fData",
                 signature = signature(
                     object = "pSet",
                     value = "data.frame"),
                 function(object, value) {
                     if (!identical(featureNames(object), rownames(value)))
                         stop("Feature names and rownames of the new fData are not identical.")
                     fd <- featureData(object)
                     pData(fd) <- value
                     object@featureData <- fd
                     if (validObject(object))
                         return(object)
                 })

setMethod("fvarMetadata",
          signature = signature(object = "pSet"),
          function(object) varMetadata(featureData(object)))
setMethod("fvarLabels",
          signature = signature(object = "pSet"),
          function(object) varLabels(featureData(object)))

setMethod("experimentData", signature(object = "pSet"),
          function(object) object@experimentData)

setMethod("msInfo", "pSet",
          function(object) msInfo(experimentData(object)))

setMethod("expinfo", "pSet",
          function(object) expinfo(experimentData(object)))

setMethod("exptitle", "pSet",
          function(object) exptitle(experimentData(object)))

setMethod("expemail", "pSet",
          function(object) expemail(experimentData(object)))

setMethod("ionSource", "pSet",
          function(object) ionSource(experimentData(object)))

setMethod("ionSourceDetails", "pSet",
          function(object) ionSourceDetails(experimentData(object)))

setMethod("analyser", "pSet",
          function(object) analyser(experimentData(object)))
setMethod("analyzer", "pSet",
          function(object) analyzer(experimentData(object)))
setMethod("analyzerDetails", "pSet",
          function(object) analyzerDetails(experimentData(object)))
setMethod("analyserDetails", "pSet",
          function(object) analyzerDetails(experimentData(object)))

setMethod("instrumentModel", "pSet",
          function(object) instrumentModel(experimentData(object)))
setMethod("instrumentManufacturer", "pSet",
          function(object) instrumentManufacturer(experimentData(object)))
setMethod("instrumentCustomisations", "pSet",
          function(object) instrumentCustomisations(experimentData(object)))

setMethod("detectorType", "pSet",
          function(object) detectorType(experimentData(object)))

setMethod("description", signature(object = "pSet"),
          function(object, ...) {
              experimentData(object)
          })
setMethod("notes", signature(object = "pSet"),
          function(object) otherInfo(experimentData(object)))

setMethod("pubMedIds", signature(object = "pSet"),
          function(object) pubMedIds(experimentData(object)))

setReplaceMethod("pubMedIds",
                 signature = signature(
                     object = "pSet",
                     value = "character"),
                 function(object, value) {
                     ed <- experimentData(object)
                     pubMedIds(ed) <- value
                     object@experimentData <- ed
                     return(object)
                 })

setMethod("abstract", "pSet",
          function(object) abstract(experimentData(object)))

setMethod("protocolData", "pSet",
          function(object) {
              tryCatch(object@protocolData,
                       error = function(x) {
                           phenoData(object)[, integer(0)]
                       })
          })

setMethod("processingData",
          signature(object = "pSet"),
          function(object) object@processingData)

setMethod("spectrapply", "pSet", function(object, FUN = NULL,
                                          BPPARAM = bpparam(), ...) {
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)
    if (is.null(FUN))
        return(spectra(object))
    bplapply(spectra(object), FUN = FUN, BPPARAM = BPPARAM, ...)
})

setMethod("$", "pSet", function(x, name) {
    eval(substitute(pData(x)$NAME_ARG, list(NAME_ARG = name)))
})
setReplaceMethod("$", "pSet", function(x, name, value) {
    pData(x)[[name]] <- value
    x
})

setReplaceMethod("pData", "pSet", function(object, value) {
    if (!is.data.frame(value))
        stop("'value' has to be a 'data.frame'")
    pData(object@phenoData) <- value
    object
})

setReplaceMethod("phenoData", "pSet", function(object, value) {
    object@phenoData <- value
    if (validObject(object))
        object
})

setMethod("isolationWindowLowerMz", "pSet", function(object)
    stop("isolationWindowLowerMz not available for ", class(object)))
setMethod("isolationWindowUpperMz", "pSet", function(object)
    stop("isolationWindowUpperMz not available for ", class(object)))

################################
## TODO: setReplaceMethods for
##  phenoData
##  pData
##  processingData - or may be individual elements for MSnProcess class
##  varMetadata
##  varLabels
##  fvarMetadata
##  fvarLabels
##  experimentData
##  description
##  notes
##  pubMedIds
##  abstract
##  protocolData

## - - - - - - - - - - - - - - - - - - - - - -
## Other TODO (based on eSet):

## setMethod("combine","pSet", ## rather a combine(MSnExp, )
