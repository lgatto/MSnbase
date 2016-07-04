## TODO: implement filter[Mz|Rt|MsLevel|File|filterAcquisitionNum] for
## MSnExp and OnDiskMSnExp classes

setMethod("filterMsLevel", "MSnExp",
          function(object, msLevel.) {
              if (!missing(msLevel.))
                  object <- object[msLevel(object) %in% msLevel.]
              object <- logging(object,
                                paste("Filter: select MS level(s)",
                                      paste(unique(msLevel.),
                                            collapse = " ")))
              object
          })
## filterMsLevel for OnDiskMSnExp:
## immediate filter: apply by subsetting the featureData.
setMethod("filterMsLevel", "OnDiskMSnExp",
          function(object, msLevel.) {
    if (!missing(msLevel.)) {
        fd <- subsetFeatureDataBy(fData(object), msLevel = msLevel.)
        fData(object) <- fd
        object <- logging(object,
                          paste0("Filter: select MS level(s) ",
                                 paste(unique(msLevel.),
                                       collapse = " "), "."))
    }
    return(object)
})

setMethod("filterMz", "MSnExp",
          function(object, mz, msLevel.) {
              ## TODO
          })

setMethod("filterRt", "MSnExp",
          function(object, rt, msLevel.) {
              ## TODO
          })

setMethod("filterFile", "MSnExp",
          function(object, file) {
              ## file can be a character or file index
              ## TODO
          })

setMethod("filterAcquisitionNum", "MSnExp",
          function(object, n) {
              ## TODO
          })

