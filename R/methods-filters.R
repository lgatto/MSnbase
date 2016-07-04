## TODO: implement filter[Mz|Rt|MsLevel|File|filterAcquisitionNum] for
## MSnExp and OnDiskMSnExp classes

setMethod("filterMsLevel", "MSnExp",
          function(object, msLevel.) {
              if (missing(msLevel.)) return(object)
              msLevel. <- as.numeric(msLevel.)
              object <- object[msLevel(object) %in% msLevel.]
              object <- logging(object,
                                paste("Filter: select MS level(s)",
                                      paste(unique(msLevel.),
                                            collapse = " ")))
              object
          })


setMethod("filterRt", "MSnExp",
          function(object, rt, msLevel.) {
              if (missing(rt)) return(object)
              if (length(rt) != 2 | !is.numeric(rt))
                  stop("'rt' must be a numeric of length 2")
              if (missing(msLevel.))
                  msLevel. <- sort(unique(msLevel(object)))
              msLevel. <- as.numeric(msLevel.)
              ## MS levels that need to be filtered out
              selms <- msLevel(object) %in% msLevel.
              ## keep spectra that match MS levels and that are within
              ## the retention time limits
              sel1 <- rtime(object) >= rt[1] & rtime(object) <= rt[2] & selms
              ## keep all spectra that don't match MS levels
              sel2 <- !selms
              object <- object[sel1 | sel2]
              fltmsg <- paste0("Filter: select retention time [",
                               paste0(rt, collapse = "-"),
                               "] and MS level(s), ",
                               paste(unique(msLevel.),
                                     collapse = " "))
              object <- logging(object, fltmsg)
              object
          })


setMethod("filterMz", "MSnExp",
          function(object, mz, msLevel.) {
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

