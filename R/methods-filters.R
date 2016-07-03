## TODO: implement filter[Mz|Rt|MsLevel|File|filterAcquisitionNum] for
## MSnExp and OnDiskMSnExp classes

setMethod("filterMsLevel", "MSnExp",
          function(object, msLevel.) {
              if (missing(msLevel.)) return(object)
              else object[msLevel(object) %in% msLevel.]
          })
