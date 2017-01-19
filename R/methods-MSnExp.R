##################################################################
## Methods for MSnExp class

setMethod("show", "MSnExp",
          function(object) {
              cat("MSn experiment data (\"", class(object), "\")\n", sep = "")
              sz <- object.size(object)
              if (length(object) == 0) {
                  show(processingData(object))
                  return(NULL)
              } else {
                  if (object@.cache$level > 0) {
                      msnMzRange <- object@.cache$rangeMz
                      rangePrecMz <- object@.cache$rangePrecursorMz
                      nPrecMz <- object@.cache$nPrecursorMz
                      uPrecMz <- object@.cache$uPrecursorMz
                      nrt <- object@.cache$nRtime
                      rtr <- object@.cache$rangeRtime
                      msLevels <- object@.cache$msLevels
                      nPrecScans <- object@.cache$nPrecursorScans
                      sz <- object@.cache$size
                  } else {
                      msLevels <- unique(msLevel(object))
                      sz <- sz +
                          sum(unlist(unname(eapply(object@assayData, object.size))))
                      msnRt <- unname(rtime(object))
                      nrt <- length(msnRt)
                      rtr <- range(msnRt)
                  }
              }
              cat("Object size in memory: ")
              cat(round(sz / (1024^2), 2), "Mb\n")
              cat("- - - Spectra data - - -\n")
              if (length(object) == 0) {
                  cat(" none\n")
              } else {
                  cat(" MS level(s):", msLevels, "\n")
                  cat(" Number of spectra:", length(object), "\n")
                  if (nrt > 0) {
                      cat(" MSn retention times:",
                          formatRt(rtr[1]), "-",
                          formatRt(rtr[2]), "minutes\n")
                  }
              }
              show(processingData(object))
              if (isOnDisk(object) &&
                  length(procs <- processingQueue(object)) > 0) {
                  cat("- - - Lazy processing queue content  - - -\n")
                  for (i in 1:length(procs))
                      cat(" o ", procs[[i]]@FUN, "\n")
              }
              cat("- - - Meta data  - - -\n")
              Biobase:::.showAnnotatedDataFrame(phenoData(object),
                                                labels=list(object="phenoData"))
              cat("Loaded from:\n")
              f <- basename(processingData(object)@files)
              nf <- length(f)
              if (nf > 0) {
                  if (nf < 3) {
                      cat(paste0("  ", f, collapse = ", "), "\n")
                  } else {
                      cat("  [1]", paste(f[1], collapse = ", "))
                      cat("...")
                      cat("  [", nf, "] ", paste(f[nf], collapse = ", "),
                          "\n", sep = "")
                      cat("  Use 'fileNames(.)' to see all files.\n")
                  }
              } else {
                  cat(" none\n")
              }
              Biobase:::.showAnnotatedDataFrame(protocolData(object),
                                                labels=list(object="protocolData"))
              Biobase:::.showAnnotatedDataFrame(featureData(object),
                                                labels=list(
                                                    object="featureData",
                                                    sampleNames="featureNames",
                                                    varLabels="fvarLabels",
                                                    varMetadata="fvarMetadata"))
              cat("experimentData: use 'experimentData(object)'\n")
              pmids <- pubMedIds(object)
              if (length(pmids) > 0 && all(pmids != ""))
                  cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
              invisible(NULL)
          })


setMethod("plot", c("MSnExp","missing"),
          function(x, y ,...) plot_MSnExp(x, ...))

setMethod("plot2d", c("MSnExp"),
          function(object, z, alpha = 1/3, plot = TRUE)
              plot2d.header(header(object), z, alpha, plot))

setMethod("plot2d", c("data.frame"),
          function(object, z, alpha = 1/3, plot = TRUE)
              plot2d.header(object, z, alpha, plot))

setMethod("plotDensity", c("MSnExp"),
          function(object, z, log = FALSE, plot = TRUE)
              plotDensity.header(header(object), z, log, plot))

setMethod("plotDensity", c("data.frame"),
          function(object, z, log = FALSE, plot = TRUE)
              plotDensity.header(object, z, log, plot))

setMethod("plotMzDelta", c("MSnExp"),
          function(object, reporters = NULL,
                   subset,
                   percentage = 0.1,
                   precMz = NULL,
                   precMzWidth = 2,
                   bw = 1,
                   xlim = c(40,200),
                   withLabels = TRUE,
                   size = 2.5,
                   plot = TRUE,
                   verbose = isMSnbaseVerbose()) {
              if (!missing(subset)) {
                  if (subset <= 0 | subset >= 1) {
                      warning('subset must be in ]0, 1[. Ignoring ',
                              subset, '.', immediate. = TRUE)
                  } else {
                      n <- length(object)
                      .subset <- sample(n, ceiling(n * subset))
                      object <- object[.subset]
                      if (verbose)
                          message("Subset to ", length(object),
                                  " spectra.")
                  }
              }
              plotMzDelta_MSnExp(object, reporters, percentage,
                                 precMz, precMzWidth,bw,
                                 xlim, withLabels, size,
                                 plot, verbose)
          })

setMethod("clean",
          signature=signature("MSnExp"),
          function(object, all = FALSE, verbose = isMSnbaseVerbose()) {
              clean_MSnExp(object, all, verbose)
          })

setMethod("removePeaks",signature("MSnExp"),
          function(object, t, verbose = isMSnbaseVerbose())
              removePeaks_MSnExp(object, t, verbose))


setMethod("trimMz",
          signature = signature("MSnExp", "numeric"),
          function(object, mzlim, msLevel.) {
              .Deprecated(new = "filterMz")
              return(filterMz(object, mz = mzlim, msLevel.))
          })

setMethod("quantify",
          signature = signature("MSnExp"),
          function(object,
                   method = c(
                       "trapezoidation", "max", "sum",
                       "SI", "SIgi", "SIn",
                       "SAF", "NSAF",
                       "count"),
                   reporters,
                   strict = FALSE,
                   parallel, ## replaced by BPPARAM
                   BPPARAM,
                   qual = TRUE,
                   verbose = isMSnbaseVerbose(),
                   ...) {
              if (!missing(parallel))
                  message("Please use BPPARAM to set a parallel framework.")
              method <- match.arg(method)
              ## this assumes that if first spectrum has msLevel > 1, all have
              if (.firstMsLevel(object) < 2)
                  stop("MS1-level quantification not implemented yet.")
              ## MS2 isobaric
              if (method %in% c("trapezoidation", "max", "sum")) {
                  if (!inherits(reporters, "ReporterIons"))
                      stop("Argument 'reporters' must inherit from 'ReporterIons' class.")
                  if (missing(BPPARAM)) {
                      BPPARAM <- bpparam()
                      if (verbose)
                          message("Using default parallel backend: ",
                                  class(BPPARAM)[1])
                  }
                  quantify_MSnExp(object, method, reporters, strict,
                                  BPPARAM, qual, verbose)
              } else if (method == "count") {
                  count_MSnSet(object)
              } else {
                  ## the following assumes that the appropriate fcols are available
                  object <- utils.removeNoIdAndMultipleAssignments(object)
                  if (method %in% c("SI", "SIgi", "SIn")) SI(object, method, ...)
                  else SAF(object, method, ...)
              }
          })

setMethod("curveStats","MSnExp",
          function(object, reporters, verbose = isMSnbaseVerbose()) {
              ifelse(verbose,progress <- "text",progress <- "none")
              l <- llply(object@spectra, curveStats, reporters, .progress = progress)
              qdfr <- l[[1]]
              for (i in 2:length(l))
                  qdfr <- rbind(qdfr,l[[i]])
              return(qdfr)
          })

setMethod("extractPrecSpectra",
          signature=signature(object="MSnExp",prec="numeric"),
          function(object,prec) extractPrecSpectra_MSnExp(object,prec))

setMethod("normalize", "MSnExp",
          function(object, method = c("max","sum"),...) {
              normalise_MSnExp(object, method = match.arg(method))
          })

normalise <- normalize

setMethod("bin", "MSnExp",
          function(object, binSize = 1, verbose = isMSnbaseVerbose()) {
              bin_MSnExp(object, binSize = binSize, verbose = verbose)
          })

setMethod("compareSpectra", c("MSnExp", "missing"),
          function(object1, fun = c("common", "cor", "dotproduct"), ...) {
              compare_MSnExp(object1, fun = match.arg(fun), ...)
          })

setMethod("pickPeaks", "MSnExp",
          function(object, halfWindowSize = 3L,
                   method = c("MAD", "SuperSmoother"),
                   SNR = 0L, ...) {
              pickPeaks_MSnExp(object, halfWindowSize = halfWindowSize,
                               method = match.arg(method), SNR = SNR, ...)
          })


setMethod("estimateNoise", "MSnExp",
          function(object, method = c("MAD", "SuperSmoother"), ...) {
              lapply(spectra(object),
                     estimateNoise_Spectrum,
                     method = match.arg(method), ...)
          })


setMethod("smooth", "MSnExp",
          function(x, method = c("SavitzkyGolay", "MovingAverage"),
                   halfWindowSize = 2L, verbose = isMSnbaseVerbose(), ...) {
              smooth_MSnExp(x, method = match.arg(method),
                            halfWindowSize = halfWindowSize, verbose = verbose,
                            ...)
          })

setMethod("removeReporters", "MSnExp",
          function(object, reporters = NULL, clean = FALSE,
                   verbose = isMSnbaseVerbose()) {
              if (is.null(reporters))
                  return(object)
              ## Throw an error if only MS1 spectra present.
              if (all(msLevel(object) == 1))
                  stop("No reporters to remove for MS1 spectra.")
              removeReporters_MSnExp(object, reporters, clean, verbose)
          })

setMethod("addIdentificationData", c("MSnExp", "character"),
          function(object, id,
                   fcol = c("spectrum.file", "acquisition.number"),
                   icol = c("spectrumFile", "acquisitionnum"),
                   verbose = isMSnbaseVerbose()) {
              addIdentificationData(object, id = mzID(id, verbose = verbose),
                                    fcol = fcol, icol = icol)
          })

setMethod("addIdentificationData", c("MSnExp", "mzIDClasses"),
          function(object, id,
                   fcol = c("spectrum.file", "acquisition.number"),
                   icol = c("spectrumFile", "acquisitionnum"), ...) {
              addIdentificationData(object, id = flatten(id),
                                    fcol = fcol, icol = icol)
          })

setMethod("addIdentificationData", c("MSnExp", "data.frame"),
          function(object, id,
                   fcol = c("spectrum.file", "acquisition.number"),
                   icol = c("spectrumFile", "acquisitionnum"), ...) {
              ## we temporaly add the spectrum.file/acquisition.number information
              ## to our fData data.frame because
              ## utils.mergeSpectraAndIdentificationData needs this information
              ## for matching (it is present in MSnSet)
              fd <- fData(object)

              if (!nrow(fd))
                  stop("No feature data found.")

              fd$spectrum.file <- basename(fileNames(object)[fromFile(object)])
              fd$acquisition.number <- acquisitionNum(object)

              fd <- utils.mergeSpectraAndIdentificationData(fd, id,
                                                            fcol = fcol,
                                                            icol = icol)

              ## after adding the identification data we remove the
              ## temporary data to avoid duplication and problems in quantify
              cn <- colnames(fd)
              keep <- cn[!(cn %in% c("spectrum.file", "acquisition.number"))]
              fData(object)[, keep] <- fd[, keep, drop=FALSE]

              if (validObject(object))
                  return(object)
          })

setMethod("removeNoId", "MSnExp",
          function(object, fcol = "pepseq", keep=NULL)
              utils.removeNoId(object, fcol, keep))

setMethod("removeMultipleAssignment", "MSnExp",
          function(object, fcol = "nprot")
              utils.removeMultipleAssignment(object, fcol))


setMethod("idSummary", "MSnExp",
          function(object) {
              ## we temporaly add the spectrumFile information
              ## to our fData data.frame because utils.idSummary
              ## needs this information for matching
              fd <- fData(object)
              fd$spectrumFile <- basename(fileNames(object)[fromFile(object)])
              return(utils.idSummary(fd))
          })

setMethod("isCentroided", "MSnExp",
          function(object, ..., verbose = isMSnbaseVerbose()) {
              pkl <- lapply(spectra(object), as.data.frame)
              ctrd <- lapply(pkl, .isCentroided, ...)
              ctrd <- unlist(ctrd, use.names = FALSE)
              if (verbose) print(table(ctrd, msLevel(object)))
              ctrd
          })

setMethod("isolationWindow", "MSnExp",
          function(object, ...) {
              ## TODO - add to fData
              mzR::isolationWindow(fileNames(object), ...)
          })
