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
              if (missing(file)) return(object)
              if (is.character(file)) {
                  file <- match(file, basename(fileNames(object)))
              }
              ## This will not work if we want to get the files in a different
              ## order (i.e. c(3, 1, 2, 5))
              file <- sort(unique(file))
              ## Sub-set spectra/featureData; all further sub-setting of
              ## processingData, experimentData and phenoData is done in [
              object <- object[fromFile(object) %in% file]
              object <- logging(object,
                                paste0("Filter: select file(s) ",
                                       paste0(file, collapse = ", "), "."))
              return(object)
          })

setMethod("filterAcquisitionNum", "MSnExp",
          function(object, n, file) {
              if (missing(n)) return(object)
              if (!is.integer(n)) stop("Argument 'n' has to be an integer",
                                       " representing the acquisition",
                                       " number(s) for sub-setting.")
              if (missing(file))
                  file <- sort(unique(fromFile(object)))
              ## Select on files.
              selFile <- fromFile(object) %in% file
              ## Select those matching the acquisition number and file.
              selAcqN <- acquisitionNum(object) %in% n & selFile
              if (!any(selAcqN))
                  warning("No spectra with the specified acquisition",
                          " number(s) found.")
              ## Subset: those from the selected files matching the acquisition
              ## num and all spectra from all other files.
              object <- object[selAcqN | !selFile]
              object <- logging(object, paste0("Filter: select by ", length(n),
                                               " acquisition number(s) in ",
                                               length(file), " file(s)."))
              return(object)
          })

