setMethod("initialize",
          signature(.Object="pSet"),
          function(.Object, ..., .cache) {
            if (missing(.cache)) {
              .cache <- new.env()
              assign("level",0,.cache)
              lockEnvironment(.cache)
            }
            callNextMethod(.Object, ..., .cache = .cache)
          })

setValidity("pSet", function(object) {
  msg <- validMsg(NULL, NULL)
  if (!all(sapply(assayData(object),function(x) inherits(x,"Spectrum"))))
    msg <- validMsg(msg,
                    "assayData must contain 'Spectrum'' objects.")
  msl <- msLevel(object)
  if (length(unique(msl))>1) 
    warning(paste("Different MS levels in ",class(object),
                  " object:",unique(msl)))
  ## checking number of spectra in assayData and
  ##          number of features in featureData
  nspectra  <- length(assayData(object)) 
  nfeatures <- nrow(featureData(object)) 
  if (nspectra != nfeatures)
    msg <- validMsg(msg,"unequal number of spectra in assayData and features in featureData.")
  if (length(spectra(object)) != length(ls(assayData(object))))
    msg <- validMsg(msg,"object size inconsistence using assayData() and spectra() methods.")  
  if (!identical(featureNames(object), ## obtained as ls(assayData(object))
                 featureNames(featureData(object))))
    msg <- validMsg(msg,"featureNames differ between assayData and featureData.")
  ## checking number of files in phenoData and
  ##          number of files in assayData
  nfilespData   <- length(processingData(object)@files)
  nfilesSpectra <- length(unique(eapply(assayData(object),fromFile)))
  if (nfilespData != nfilesSpectra)
    msg <- validMsg(msg, "unequal number of files in assayData and phenoData.")  
  ## protocolData not checked yet - depends very much
  ## on type of assay (MS1, MS2 quant, reporter ions, ...)
  if (!cacheEnvIsLocked(object))
    msg <- validMsg(msg,"'.cache' environment is not locked.")
  if (!exists("level",envir=object@.cache))
    msg <- validMsg(msg,"'.cache' level not defined.")
  if (is.null(msg)) TRUE else msg
})


setMethod("[","pSet",
          function(x,i,j="missing",drop="missing") {
            if (!(is.logical(i) | is.numeric(i)))
              stop("subsetting works only with numeric or logical")
            if (is.numeric(i)) {
              if (max(i)>length(x) | min(i)<1)
                stop("subscript out of bounds")
            }
            whichElements <- ls(assayData(x))[i]
            x@assayData <- list2env(mget(whichElements,assayData(x)))
            x@featureData <- featureData(x)[i,]
            if (is.logical(i)) {
              x@processingData@processing <-
                c(processingData(x)@processing,
                  paste("Data [logically] subsetted ",sum(i)," spectra: ",date(),sep=""))
            } else if (is.numeric(i)) {
              x@processingData@processing <-
                c(processingData(x)@processing,
                  paste("Data [numerically] subsetted ",length(i)," spectra: ",date(),sep=""))
            } else {
              x@processingData@processing <-
                c(processingData(x)@processing,
                  paste("Data subsetted ",i,": ",date(),sep=""))
            }
            if (x@.cache$level > 0)
              x@.cache <- setCacheEnv(assayData(x), x@.cache$level)
            if (validObject(x))
              return(x)
          })


setMethod("[[","pSet",
          function(x,i,j="missing",drop="missing") {
            if (length(i)!=1)
              stop("subscript out of bounds")
            if (!is.character(i)) 
              i <- featureNames(x)[i]
            return(get(i,envir=assayData(x)))
          })

setMethod("precursorMz","pSet",
          function(object) {
            ## this assumes that if first spectrum
            ## has msLevel>1, all have
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), precursorMz))
            stop("No precursor MZ value for MS1 spectra.")
          })

setMethod("precScanNum","pSet",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(unlist(sapply(spectra(object), precScanNum)))
            stop("This experiment contains MS1 spectra.")
          })


setMethod("tic","pSet",
          function(object) sapply(spectra(object),tic))

setMethod("precursorCharge","pSet",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), precursorCharge))
            stop("No precursor MZ value for MS1 spectra.")
          })

setMethod("precursorIntensity","pSet",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object), precursorIntensity))
            stop("No precursor data for MS1 spectra.")
          })

setMethod("acquisitionNum","pSet",
          function(object) sapply(spectra(object), acquisitionNum))

setMethod("rtime","pSet",
          function(object) sapply(spectra(object),rtime))

setMethod("centroided","pSet",
          function(object) sapply(spectra(object),centroided))

setReplaceMethod("centroided",
                 signature(object="pSet",
                           value="logical"),
                 function(object, value) {
                   if (length(object) != length(value))
                     stop("Length of replacement value is different than number of spectra.")
                   sl <- spectra(object)
                   for (i in 1:length(sl))
                     centroided(sl[[i]]) <- value[i]
                   object@assayData <- as.environment(sl)
                   if (validObject(object))
                     return(object)
                 })

setMethod("peaksCount",
          signature("pSet","missing"),
          function(object) sapply(spectra(object),peaksCount))

setMethod("peaksCount",
          signature("pSet","numeric"),          
          function(object, scans) {
            if (length(scans)==1)
              return(peaksCount(object[[scans]]))           
            sapply(spectra(object)[scans],peaksCount)
          })

setMethod("msLevel","pSet",
          function(object) sapply(spectra(object),msLevel))

setMethod("collisionEnergy","pSet",
          function(object) {
            if (msLevel(object)[1]>1) 
              return(sapply(spectra(object),collisionEnergy))
            stop("No collision energy for MS1 spectra.")
          })

setMethod("intensity","pSet",
          function(object) lapply(spectra(object),intensity))

setMethod("mz","pSet",
          function(object) lapply(spectra(object),mz))

setMethod("polarity","pSet",
          function(object) {
            if (msLevel(object)[1]==1) 
              return(sapply(spectra(object), polarity))
            stop("No polarity for MS2 spectra.")
          })

setMethod("fromFile","pSet",
          function(object) return(sapply(spectra(object),fromFile)))

setMethod("header",
          signature("pSet","missing"),
          function(object) {
            if (any(msLevel(object) < 2))
              stop("header() only works for MS levels > 1.")
            tbl <- table(fromFile(object))
            idx <- as.numeric(unlist(apply(tbl, 1, function(x) 1:x)))
            return(data.frame(cbind(index = idx,
                                    file = fromFile(object),
                                    retention.time = rtime(object),
                                    precursor.mz = precursorMz(object),
                                    precursor.intensity = precursorIntensity(object),
                                    charge = precursorCharge(object),
                                    peaks.count = peaksCount(object),
                                    tic = tic(object),
                                    ms.level = msLevel(object),
                                    acquisition.number = acquisitionNum(object),
                                    collision.energy = collisionEnergy(object))))
          })

setMethod("header",
          signature=c("pSet","numeric"),
          function(object, scans) {
            if (any(msLevel(object) < 2))
              stop("header() only works for MS levels > 1.")
            tbl <- table(fromFile(object))
            idx <- as.numeric(unlist(apply(tbl, 1, function(x) 1:x)))
            ## OK for length(scans) > 1 -- slow for 1 scan
            hdfr <- data.frame(cbind(index = idx,
                                     file = fromFile(object),
                                     retention.time = rtime(object),
                                     precursor.mz = precursorMz(object),
                                     peaks.count = peaksCount(object),
                                     tic = tic(object),
                                     ms.level = msLevel(object),
                                     charge = precursorCharge(object),
                                     collision.energy = collisionEnergy(object)))          
            return(hdfr[scans, ])
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

setMethod("spectra","MSnExp",function(object) {
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

## No proteomicsData slot anymore since MSnbase 0.2.0 - experimentData is not MIAPE
## setMethod("proteomicsData","MSnExp",function(object) object@proteomicsData)
## setMethod("proteomicsData<-","MSnExp",
##           function(object,value="MIAPE") object@proteomicsData <- value)
## setReplaceMethod("proteomicsData",
##                  signature(object="MSnExp",
##                            value="MIAPE"),
##                  function(object, value) {
##                    object@proteomicsData = value
##                    if (validObject(object))
##                      return(object)
##                  })

setMethod("sampleNames",
          signature(object="pSet"),
          function(object) sampleNames(phenoData(object)))

setMethod("fileNames",
          signature(object="pSet"),
          function(object) processingData(object)@files)

## setReplaceMethod("sampleNames",
##                  signature=signature(object="pSet", value="character"),
##                  function(object, value) {
##                      pd <- phenoData(object)
##                      sampleNames(pd) <- value
##                      prd <- protocolData(object)
##                      if (nrow(prd) == 0) {
##                        prd <- pd[,integer(0)]
##                      } else {
##                        sampleNames(prd) <- value
##                      }
##                      object@phenoData <- pd
##                      object@protocolData <- prd
##                      if (validObject(object))
##                        return(object)
##                  })

setMethod("featureNames",
          signature=signature(object="pSet"),
          function(object) ls(assayData(object)))

setMethod("phenoData", "pSet", function(object) object@phenoData)
setMethod("pData", "pSet", function(object) pData(phenoData(object)))
setMethod("varMetadata",
          signature=signature(object="pSet"),
          function(object) varMetadata(phenoData(object)))
setMethod("varLabels",
          signature=signature(object="pSet"),
          function(object) varLabels(phenoData(object)))
setMethod("featureData",
          signature(object="pSet"),
          function(object) object@featureData)

setReplaceMethod("featureData",
                 signature=signature(
                   object="pSet",
                   value="AnnotatedDataFrame"),
                 function(object, value) {
                   object@featureData <- value
                   if (validObject(object))
                     return(object)
                 })

setMethod("fData",
          signature=signature(object="pSet"),
          function(object) pData(featureData(object)))

setReplaceMethod("fData",
                 signature=signature(
                   object="pSet",
                   value="data.frame"),
                 function(object, value) {
                     fd <- featureData(object)
                     pData(fd) <- value
                     object@featureData <- fd
                     if (validObject(object))
                       return(object)
                 })

setMethod("fvarMetadata",
          signature=signature(object="pSet"),
          function(object) varMetadata(featureData(object)))
setMethod("fvarLabels",
          signature=signature(object="pSet"),
          function(object) varLabels(featureData(object)))
setMethod("experimentData", signature(object="pSet"), function(object) object@experimentData)

setMethod("msInfo","pSet",
          function(object) msInfo(experimentData(object)))

setMethod("description", signature(object="pSet"),
          function(object, ...) {
            experimentData(object)
          })
setMethod("notes", signature(object="pSet"),
          function(object) otherInfo(experimentData(object)))

setMethod("pubMedIds", signature(object="pSet"),
          function(object) pubMedIds(experimentData(object)))

setReplaceMethod("pubMedIds",
                 signature=signature(
                   object="pSet",
                   value="character"),
                 function(object, value) {
                     ed <- experimentData(object)
                     pubMedIds(ed) <- value
                     object@experimentData <- ed
                     return(object)
                   })


setMethod("abstract", "pSet", function(object) abstract(experimentData(object)))

setMethod("protocolData", "pSet",
          function(object) {
            tryCatch(object@protocolData,
                     error = function(x) {
                       phenoData(object)[,integer(0)]
                     })
          })
setMethod("processingData",
          signature(object="pSet"),
          function(object) object@processingData)


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

## setMethod("combine","pSet", ## rather a combine(MSnSet, )

## setReplaceMethod("featureNames",
##                  signature=signature(object="pSet", value="ANY"),
##                  function(object, value) {
##                    fd <- featureData(object)
##                    featureNames(fd) <- value
##                    ad <- assayData(object)
##                    featureNames(ad) <- value
##                    object@featureData <- fd
##                    return(object)
##                  })

