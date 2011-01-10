## No initialize method for pSet -- use constructor

setValidity("pSet", function(object) {
  msg <- validMsg(NULL, NULL)
  ## checking number of spectra in assayData and
  ##          number of features in featureData
  nspectra  <- length(assayData(object)) 
  nfeatures <- nrow(featureData(object)) 
  if (nspectra != nfeatures)
    msg <- validMsg(msg, "unequal number of spectra in assayData and features in featureData")
  ## checking number of files in phenoData and
  ##          number of files in assayData
  nfilespData   <- length(object@process@files)
  nfilesSpectra <- length(unique(eapply(assayData(object),fromFile)))
  if (nfilespData != nfilesSpectra)
    msg <- validMsg(msg, "unequal number of files in assayData and phenoData")  
  ## protocolData not checked yet - depends very much
  ## on type of assay (MS1, MS2 quant, reporter ions, ...)
  if (is.null(msg)) TRUE else msg
})

## returns the dimensions of PhenoData 
setMethod("dim", "pSet", function(x) dim(pData(x)))

## returns the number of spectra in the AssayData env
setMethod("length", "pSet", function(x) length(assayData(x)))

setMethod("assayData", "pSet", function(object) object@assayData)

## setReplaceMethod("assayData",
##                  signature=signature(
##                    object="pSet",
##                    value="AssayData"),
##                  function(object, value) {
##                    object@assayData <- value ## may be lock env?                   
##                    return(object)
##                  })

setMethod("sampleNames",
          signature(object="pSet"),
          function(object) sampleNames(phenoData(object)))

setMethod("fileNames",
          signature(object="pSet"),
          function(object) sampleNames(object))

setReplaceMethod("sampleNames",
                 signature=signature(object="pSet", value="character"),
                 function(object, value) {
                     pd <- phenoData(object)
                     sampleNames(pd) <- value
                     prd <- protocolData(object)
                     if (nrow(prd) == 0) {
                       prd <- pd[,integer(0)]
                     } else {
                       sampleNames(prd) <- value
                     }
                     object@phenoData <- pd
                     object@protocolData <- prd
                     if (validObject(object))
                       return(object)
                 })

setReplaceMethod("fileNames",
                 signature=signature(object="pSet", value="character"),
                 function(object, value) sampleNames(object) <- value)


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
setMethod("fData",
          signature=signature(object="pSet"),
          function(object) pData(featureData(object)))
setMethod("fvarMetadata",
          signature=signature(object="pSet"),
          function(object) varMetadata(featureData(object)))
setMethod("fvarLabels",
          signature=signature(object="pSet"),
          function(object) varLabels(featureData(object)))
setMethod("experimentData", signature(object="pSet"), function(object) object@experimentData)
setMethod("description", signature(object="pSet"),
          function(object, ...) {
            experimentData(object)
          })
setMethod("notes", signature(object="pSet"),
          function(object) otherInfo(experimentData(object)))
setMethod("pubMedIds", signature(object="pSet"),
          function(object) pubMedIds(experimentData(object)))
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
          function(object) object@process)

################################
## TODO: setReplaceMethods for
##  phenoData
##  pData
##  processingData - or may be individual elements for MSnProcess class
##  varMetadata
##  varLabels
##  featureData
##  fData
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

## setMethod("combine","pSet",

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


## setMethod("[", "pSet",
## setMethod("$", "pSet",
## setReplaceMethod("$", "pSet",
## setMethod("[[", "pSet", 
## setReplaceMethod("[[", "pSet",

