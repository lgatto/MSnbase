############################################################
## Methods for OnDiskMSnExp objects.

setMethod("initialize",
          signature(.Object="OnDiskMSnExp"),
          function(.Object, ...){
              ## if(!any(names(list(...)) == "onDisk")){
              ##     .Object@onDisk=TRUE
              ## }
              callNextMethod()
          })

############################################################
## msLevel
##
## Extract the msLevel info for all spectra in an OnDiskMSnExp
## object. In contrast to the MSnExp we're not getting that
## from the individual Spectrum objects in @assayData, but
## from the featureData.
setMethod("msLevel", "OnDiskMSnExp", function(object){
    msl <- featureData(object)$msLevel
    ## That should not happen!
    if(is.null(msl))
        stop("The 'OnDiskMSnExp' does not have msLevel information in the featureData!")
    names(msl) <- featureNames(object)
    return(msl)
})

############################################################
## fromFile
##
## Extract the index from which file the spectra in the data set
## derive from.
setMethod("fromFile", "OnDiskMSnExp", function(object){
    fidx <- fData(object)$fileIdx
    names(fidx) <- featureNames(object)
    return(fidx)
})

