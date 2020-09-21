## TODO: implement filter[Mz|Rt|MsLevel|File|filterAcquisitionNum] for
## MSnExp and OnDiskMSnExp classes

setMethod("filterMsLevel", "MSnExp",
          function(object, msLevel.) {
              if (missing(msLevel.)) return(object)
              msLevel. <- as.numeric(msLevel.)
              object <- object[msLevel(object) %in% msLevel.]
              object <- logging(object,
                                paste0("Filter: select MS level(s) ",
                                      paste(unique(msLevel.),
                                            collapse = " "), "."))
              object
          })

setMethod("filterPolarity", "MSnExp",
          function(object, polarity.) {
              if (missing(polarity.)) return(object)
              msLevel. <- as.numeric(polarity.)
              object <- object[polarity(object) %in% polarity.]
              object <- logging(object,
                                paste0("Filter: select polarity ",
                                polarity., "."))
              object
          })


setMethod("filterRt", "MSnExp",
          function(object, rt, msLevel.) {
              if (missing(rt)) return(object)
              if (length(rt) != 2 | !is.numeric(rt))
                  stop("'rt' must be a numeric of length 2")
              if (missing(msLevel.))
                  msLevel. <- base::sort(unique(msLevel(object)))
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

############################################################
## filterMz for MSnExp
setMethod("filterMz", "MSnExp",
          function(object, mz, msLevel., ...) {
              if (missing(mz))
                  return(object)
              if (!is.numeric(mz) & length(mz) != 2)
                  stop("'mz' must be a numeric of length 2!")
              if (missing(msLevel.)) {
                  msLevel. <- base::sort(unique(msLevel(object)))
              } else {
                  if (!is.numeric(msLevel.))
                      stop("'msLevel' must be numeric!")
              }
              ## Note: the msLevel. argument is passed down to the
              ## filterMz_Spectrum function.
              filtered <- eapply(assayData(object), filterMz, mz = mz,
                                 msLevel. = msLevel., ...)
              object@assayData <- list2env(filtered)
              trmd <- object@processingData@trimmed
              ifelse(length(trmd) == 0,
                     object@processingData@trimmed <- mz,
                     object@processingData@trimmed <- c(max(trmd[1], mz[1]),
                                                        min(trmd[2], mz[2])))
              object@processingData@processing <- c(object@processingData@processing,
                                                    paste0("Filter: trim MZ [",
                                                          object@processingData@trimmed[1],
                                                          "..",object@processingData@trimmed[2],
                                                          "] on MS level(s) ",
                                                          paste(msLevel., sep = ", "), "."))
              if (object@.cache$level > 0) {
                  hd <- header(object)
                  hd$peaks.count <- peaksCount(object)
                  hd$ionCount <- ionCount(object)
                  object@.cache <- setCacheEnv(list(assaydata = assayData(object),
                                                    hd = hd),
                                               object@.cache$level)
              }
              if (validObject(object))
                  return(object)
          })
## filterMz for OnDiskMSnExp
setMethod("filterMz", "OnDiskMSnExp",
          function(object, mz, msLevel., ...) {
              if (missing(mz))
                  return(object)
              if (!is.numeric(mz) & length(mz) != 2)
                  stop("'mz' must be a numeric of length 2!")
              if (missing(msLevel.)) {
                  msLevel. <- base::sort(unique(msLevel(object)))
              } else {
                  if (!is.numeric(msLevel.))
                      stop("'msLevel' must be numeric!")
              }
              ps <- ProcessingStep("filterMz", list(mz = mz,
                                                    msLevel. = msLevel., ...))
              object@spectraProcessingQueue <- c(object@spectraProcessingQueue,
                                                 list(ps))
              trmd <- object@processingData@trimmed
              ifelse(length(trmd)==0,
                     object@processingData@trimmed <- mz,
                     object@processingData@trimmed <- c(max(trmd[1],mz[1]),
                                                        min(trmd[2],mz[2])))
              object@processingData@processing <- c(object@processingData@processing,
                                                    paste0("Filter: trim MZ [",
                                                          object@processingData@trimmed[1],
                                                          "..",object@processingData@trimmed[2],
                                                          "] on MS level(s) ",
                                                          paste(msLevel., sep = ", "), "."))
              return(object)
          })


setMethod("filterFile", "MSnExp",
          function(object, file) {
              ## file can be a character or file index
              if (missing(file)) return(object)
              if (is.character(file)) {
                  file <- base::match(file, basename(fileNames(object)))
              }
              ## This will not work if we want to get the files in a different
              ## order (i.e. c(3, 1, 2, 5))
              file <- base::sort(unique(file))
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
                  file <- base::sort(unique(fromFile(object)))
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


setMethod("filterEmptySpectra", "MSnExp",
          function(object, ...) {
              emptyspecs <- unlist(spectrapply(object, isEmpty))
              object <- object[!emptyspecs, ]
              msg <- paste("Removed", sum(emptyspecs), "empty spectra.")
              object <- logging(object, msg)
              object
          })

setMethod("filterPrecursorScan", "MSnExp",
          function(object, acquisitionNum, ...) {
            object <- object[.filterSpectraHierarchy(fData(object),
                                                     acquisitionNum)]
            object <- logging(
                object,
                paste(
                    "Filter: select parent/children scans for",
                    paste0(acquisitionNum, collapse=", ")
                )
            )
            object
        })

#' @description
#'
#' Filter `object` reducing to MS data (spectra) for which the precursor m/z
#' matches the provided m/z accepting a small difference which can be defined
#' with the `ppm` parameter.
#'
#' @param object `MSnExp` object.
#'
#' @param mz `numeric(1)` with the m/z value to filter the object.
#'
#' @param ppm `numeric(1)` defining the accepted difference between the
#'     provided m/z and the spectrum's m/z in parts per million.
#'
#' @md
#'
#' @noRd
#'
#' @author Johannes Rainer
setMethod("filterPrecursorMz", "MSnExp", function(object, mz, ppm = 10) {
    if (missing(mz))
        return(object)
    if (length(mz) > 1)
        stop("'mz' is expected to be a single m/z value", call. = FALSE)
    mz_ppm <- ppm * mz / 1e6
    keep <- which(precursorMz(object) >= (mz - mz_ppm) &
                 precursorMz(object) <= (mz + mz_ppm))
    object <- object[keep]
    msg <- paste0("Filter: select by precursor m/z: ", mz, ".")
    object <- logging(object, msg)
    object
})

#' @description
#'
#' Filter the object based on the isolation window. The resulting object will
#' contain all spectra with isolation windows containing the specified m/z.
#'
#' @param object `MSnExp` object.
#'
#' @param mz `numeric(1)` with the m/z tha should be within the isolation
#'     window's m/z range.
#'
#' @md
#'
#' @noRd
#'
#' @author Johannes Rainer
setMethod("filterIsolationWindow", "MSnExp", function(object, mz, ...){
    if (missing(mz))
        return(object)
    if (length(mz) > 1)
        stop("'mz' is expected to be a single m/z value", call. = FALSE)
    keep <- which(isolationWindowLowerMz(object) <= mz &
                  isolationWindowUpperMz(object) >= mz)
    object <- object[keep]
    msg <- paste0("Filter isolation windows containing : ", mz, ".")
    object <- logging(object, msg)
    object
})
