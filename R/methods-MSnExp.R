##################################################################
## Methods for MSnExp class

setMethod("show",
          signature=signature(object="MSnExp"),
          function(object) {
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
              if (all(msLevel(object) > 1)) {
                msnPrecMz <- unname(eapply(assayData(object),precursorMz))
                nPrecMz <- length(msnPrecMz)
                uPrecMz <- length(unique(msnPrecMz))
                rangePrecMz <- range(msnPrecMz)
                msnMzRange <- range(unname(mz(object)))
                nPrecScans <- length(unique(eapply(assayData(object),precScanNum)))
              }
              msLevels <- unique(unlist(eapply(assayData(object),msLevel)))
              sz <- sum(unlist(unname(eapply(assayData(object),object.size))))
              msnRt <- unname(eapply(assayData(object),rtime))
              nrt <- length(msnRt)
              rtr <- range(msnRt)
            }
            cat("Object of class \"",class(object),"\"\n",sep="")
            cat(" Object size in memory: ")
            if (length(assayData(object)) == 0) {
              sz <- object.size(object)
            } else {
              sz <- sz + object.size(object)
            }
            cat(round(sz/(1024^2),2),"Mb\n")
            cat("- - - Spectra data - - -\n")
            if (length(assayData(object)) == 0) {
              cat(" none\n")
            } else {
              cat(" MS level(s):",msLevels,"\n")
              if (all(msLevel(object) > 1)) {
                cat(" Number of MS1 acquisitions:",nPrecScans,"\n")
                cat(" Number of MSn scans:",length(ls(assayData(object))),"\n")
                cat(" Number of precursor ions:",nPrecMz,"\n")
                if (nPrecMz > 0) {
                  cat("",uPrecMz,"unique MZs\n")
                  cat(" Precursor MZ's:",paste(signif(rangePrecMz,5),collapse=" - "),"\n")
                }
                cat(" MSn M/Z range:",round(msnMzRange,2),"\n")
              } else {
                cat(" Number of MS1 scans:",length(spectra(object)),"\n")
              }
              if (nrt > 0) {
                cat(" MSn retention times:",formatRt(rtr[1]),"-",formatRt(rtr[2]),"minutes\n")
              }
            }
            show(processingData(object))
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


setMethod("plot",c("MSnExp","missing"),
          function(x,y,...) plot_MSnExp(x,...))

setMethod("plot2d",c("MSnExp"),
          function(object,z,alpha=1/3,plot=TRUE)
          plot2d.header(header(object),z,alpha,plot))

setMethod("plot2d",c("data.frame"),
          function(object,z,alpha=1/3,plot=TRUE)
          plot2d.header(object,z,alpha,plot))

setMethod("plotDensity",c("MSnExp"),
          function(object,z,log=FALSE,plot=TRUE)
          plotDensity.header(header(object),z,log,plot))

setMethod("plotDensity",c("data.frame"),
          function(object,z,log=FALSE,plot=TRUE)
          plotDensity.header(object,z,log,plot))

setMethod("plotMzDelta",c("MSnExp"),
          function(object, reporters=NULL,
                   subset,
                   percentage=0.1,
                   precMz=NULL,
                   precMzWidth=2,
                   bw=1,
                   xlim=c(40,200),
                   withLabels=TRUE,
                   size=2.5,
                   plot=TRUE,
                   verbose=TRUE) {
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
          function(object, all = FALSE, verbose = TRUE) clean_MSnExp(object, all, verbose))

setMethod("removePeaks",signature("MSnExp"),
          function(object, t, verbose = TRUE) removePeaks_MSnExp(object, t, verbose))


setMethod("trimMz",
          signature=signature("MSnExp","numeric"),
          function(object, mzlim, ...) {
            trimmed <- eapply(assayData(object), trimMz, mzlim, ...)
            object@assayData <- list2env(trimmed)
            trmd <- object@processingData@trimmed
            ifelse(length(trmd)==0,
                   object@processingData@trimmed <- mzlim,
                   object@processingData@trimmed <- c(max(trmd[1],mzlim[1]),
                                                      min(trmd[2],mzlim[2])))
            object@processingData@processing <- c(object@processingData@processing,
                                                  paste("MZ trimmed [",object@processingData@trimmed[1],
                                                        "..",object@processingData@trimmed[2],"]",sep=""))
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
                   verbose = TRUE,
                   ...) {
              if (!missing(parallel))
                  message("Please use BPPARAM to set a parallel framework.")
              method <- match.arg(method)
              ## this assumes that if first spectrum has msLevel > 1, all have
              if (msLevel(object)[1] < 2)
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
          function(object, reporters, verbose = TRUE) {
            ifelse(verbose,progress <- "text",progress <- "none")
            l <- llply(object@spectra, curveStats, reporters, .progress=progress)
            qdfr <- l[[1]]
            for (i in 2:length(l))
              qdfr <- rbind(qdfr,l[[i]])
            return(qdfr)
          })

setMethod("extractPrecSpectra",
          signature=signature(object="MSnExp",prec="numeric"),
          function(object,prec) extractPrecSpectra_MSnExp(object,prec))

setMethod("extractSpectra",
          signature=signature(object="MSnExp",selected="logical"),
          function(object,selected) {
            msg <- c("The 'extractSpectra' function is defunct\n",
                     "Please use the '[' subsetting operator instead.")
            .Defunct(msg=msg)
            ## extractSpectra.MSnExp(object,selected)
          })

setMethod("normalize", "MSnExp",
          function(object, method = c("max","sum"),...) {
            normalise_MSnExp(object, method = match.arg(method))
        })

normalise <- normalize

setMethod("bin", "MSnExp",
          function(object, binSize = 1, verbose = TRUE) {
            bin_MSnExp(object, binSize = binSize, verbose = verbose)
        })

setMethod("compareSpectra", c("MSnExp", "missing"),
          function(object1, fun = c("common", "cor", "dotproduct"), ...) {
            compare_MSnExp(object1, fun=match.arg(fun), ...)
        })

setMethod("pickPeaks", "MSnExp",
          function(object, halfWindowSize = 3L,
                   method = c("MAD", "SuperSmoother"),
                   SNR = 0L, ...) {
            pickPeaks_MSnExp(object, halfWindowSize = halfWindowSize,
                             method = match.arg(method), SNR = SNR, ...)
        })

setMethod("smooth", "MSnExp",
          function(x, method = c("SavitzkyGolay", "MovingAverage"),
                   halfWindowSize = 2L, verbose = TRUE, ...) {
            smooth_MSnExp(x, method = match.arg(method),
                          halfWindowSize = halfWindowSize, verbose = verbose,
                          ...)
        })

setMethod("removeReporters","MSnExp",
          function(object, reporters=NULL, clean=FALSE, verbose=TRUE) {
            if (is.null(reporters))
              return(object)
            removeReporters_MSnExp(object, reporters, clean, verbose)
        })

setMethod("addIdentificationData", "MSnExp",
          function(object, filenames, verbose = TRUE) {
            ## we temporaly add the file/acquisition.number information
            ## to our fData data.frame because utils.addIdentificationData
            ## needs this information for matching (it is present in MSnSet)
            fData(object)$file <- fromFile(object)
            fData(object)$acquisition.number <- acquisitionNum(object)
            object <- utils.addIdentificationData(object, filenames,
                                                  verbose = verbose)
            ## after adding the identification data we remove the
            ## temporary data to avoid duplication and problems in quantify
            cn <- colnames(fData(object))
            keep <- !(cn %in% c("file", "acquisition.number"))
            fData(object) <- fData(object)[, keep, drop=FALSE]
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
            ## we temporaly add the file information
            ## to our fData data.frame because utils.idSummary
            ## needs this information for matching (it is present in MSnSet)
            fd <- fData(object)
            fd$file <- fromFile(object)
            return(utils.idSummary(fd))
        })
