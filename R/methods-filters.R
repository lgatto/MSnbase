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

