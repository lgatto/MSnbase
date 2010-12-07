## No initialize method for pSet -- use constructor

## returns the dimensions of PhenoData 
setMethod("dim", "pSet", function(x) dim(pData(x)))

## returns the number of spectra in the AssayData env
setMethod("length", "pSet", function(x) length(ls(assayData(x))))

setValidity("eSet", function(object) {
  msg <- validMsg(NULL, NULL)
  dims <- dims(object)
  if (ncol(dims) > 0) {
    ## assayData
    ## featureData
    nspectra  <- length(assayData(object)) ## number of spectra in assayData
    nfeatures <- nrow(featureData(object)) ## number of features in featureData
    if (nspectra != nfeatures)
      msg <- validMsg(msg, "unequal number of spectra in assayData and features in featureData")
    ## phenoData
    nfilespData   <- nrow(pData(object))                                ## number of files in phenoData
    nfilesSpectra <- length(unique(eapply(assayData(object),fromFile))) ## number of files in assayData
    if (nfilespData != nfilesSpectra)
      msg <- validMsg(msg, "unequal number of files in assayData and phenoData")
    ## protocolData -- as in eSet
    if (dim(phenoData(object)) != dim(protocolData(object)))
      msg <- validMsg(msg, "sample numbers differ between phenoData and protocolData")
    if (!identical(sampleNames(phenoData(object)), sampleNames(protocolData(object))))
      msg <- validMsg(msg, "sampleNames differ between phenoData and protocolData")
  }
  if (is.null(msg)) TRUE else msg
})

setMethod("assayData", "pSet", function(object) object@assayData)

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

## TODO: setReplaceMethods for
##  phenoData
##  pData
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

## setReplaceMethod("assayData",
##                  signature=signature(
##                    object="pSet",
##                    value="AssayData"),
##                  function(object, value) {
##                    object@assayData <- value ## may be lock env?                   
##                    return(object)
##                  })
