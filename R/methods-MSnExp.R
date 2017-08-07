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
                   pepseq = "sequence",
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
                  count_MSnSet(object, pepseq)
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
             icol = c("spectrumFile", "acquisitionNum"),
             acc = "DatabaseAccess",
             desc = "DatabaseDescription",
             pepseq = "sequence",
             key = "spectrumID",
             decoy = "isDecoy",
             rank = "rank",
             accession = acc,
             verbose = isMSnbaseVerbose(),
             ...)
        .addCharacterIdentificationData(object, id, fcol, icol, acc,
                                        desc, pepseq, key, decoy,
                                        rank, accession, verbose, ...))

setMethod("addIdentificationData", c("MSnExp", "mzRident"),
        function(object, id,
                 fcol = c("spectrum.file", "acquisition.number"),
                 icol = c("spectrumFile", "acquisitionNum"),
                 acc = "DatabaseAccess",
                 desc = "DatabaseDescription",
                 pepseq = "sequence",
                 key = "spectrumID",
                 decoy = "isDecoy",
                 rank = "rank",
                 accession = acc,
                 verbose = isMSnbaseVerbose(),
                 ...)
            .addMzRidentIdentificationData(object, id, fcol, icol,
                                           acc, desc, pepseq, key,
                                           decoy, rank, accession,
                                           verbose, ...))

setMethod("addIdentificationData", c("MSnExp", "mzIDClasses"),
        function(object, id,
                 fcol = c("spectrum.file", "acquisition.number"),
                 icol = c("spectrumFile", "acquisitionnum"),
                 acc = "accession",
                 desc = "description",
                 pepseq = "pepseq",
                 key = "spectrumid",
                 decoy = "isdecoy",
                 rank = "rank",
                 accession = acc,
                 verbose = isMSnbaseVerbose(),
                 ...)
            .addMzIDIdentificationData(object, id, fcol, icol, acc,
                                       desc, pepseq, key, decoy, rank,
                                       accession, verbose,...))

setMethod("addIdentificationData", c("MSnExp", "data.frame"),
          function(object, id,
                   fcol = c("spectrum.file", "acquisition.number"),
                   icol, acc, desc, pepseq, key, decoy, rank,
                   accession = acc, verbose = isMSnbaseVerbose(), ...)
              .addDataFrameIdentificationData(object, id, fcol, icol,
                                              acc, desc, pepseq, key,
                                              decoy, rank, accession,
                                              verbose, ...))

setMethod("removeNoId", "MSnExp",
          function(object, fcol = "sequence", keep=NULL)
              utils.removeNoId(object, fcol, keep))

setMethod("removeMultipleAssignment", "MSnExp",
          function(object, nprot = "nprot")
              utils.removeMultipleAssignment(object, nprot))


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

setMethod("splitByFile", c("MSnExp", "factor"), function(object, f) {
    if (length(f) != length(fileNames(object)))
        stop("length of 'f' has to match the length of samples/files in 'object'.")
    idxs <- lapply(levels(f), function(z) which(f == z))
    ## Use filterFile to split them.
    res <- lapply(idxs, function(z) {
        return(filterFile(object, file = z))
    })
    names(res) <- levels(f)
    return(res)
})

#' @title Extract chromatogram object(s)
#'
#' @aliases chromatogram
#'
#' @description The \code{chromatogram} method extracts chromatogram(s) from an
#'     \code{\linkS4class{MSnExp}} or \code{\linkS4class{OnDiskMSnExp}} object.
#'     Depending on the provided parameters this can be a total ion chromatogram
#'     (TIC), a base peak chromatogram (BPC) or an extracted ion chromatogram
#'     (XIC) extracted from each sample/file.
#'
#' @details Arguments \code{rt} and \code{mz} allow to specify the MS
#'     data slice from which the chromatogram should be extracted.
#'     The parameter \code{aggregationSum} allows to specify the function to be
#'     used to aggregate the intensities across the mz range for the same
#'     retention time. Setting \code{aggregationFun = "sum"} would e.g. allow
#'     to calculate the \emph{total ion chromatogram} (TIC),
#'     \code{aggregationFun = "max"} the \emph{base peak chromatogram} (BPC).
#'     The length of the extracted \code{\link[MSnbase]{Chromatogram}} object,
#'     i.e. the number of available data points, corresponds to the number of
#'     scans/spectra measured in the specified retention time range. If in a
#'     specific scan (for a give retention time) no signal was measured in the
#'     specified mz range, a \code{NA_real_} is reported as intensity for the
#'     retention time (see Notes for more information). This can be changed
#'     using the \code{missing} parameter.
#'
#'     By default or if \code{mz} and/or \code{rt} are numeric vectors, the
#'     function extracts one \code{\link{Chromatogram}} object for each file
#'     in the \code{\linkS4class{MSnExp}} or \code{\linkS4class{OnDiskMSnExp}}
#'     object. Providing a numeric matrix with argument \code{mz} or \code{rt}
#'     enables to extract multiple chromatograms per file, one for each row in
#'     the matrix. If the number of columns of \code{mz} or \code{rt} are not
#'     equal to 2, \code{range} is called on each row of the matrix.
#' 
#' @param object For \code{chromatogram}: a \code{\linkS4class{MSnExp}} or
#'     \code{\linkS4class{OnDiskMSnExp}} object from which the chromatogram
#'     should be extracted.
#'
#' @param rt A \code{numeric(2)} or two-column \code{matrix} defining the lower
#'     and upper boundary for the retention time range/window(s) for the
#'     chromatogram(s). If a \code{matrix} is provided, a chromatogram is
#'     extracted for each row. If not specified, a chromatogram representing the
#'     full retention time range is extracted. See examples below for details.
#'
#' @param mz A \code{numeric(2)} or two-column \code{matrix} defining the
#'     mass-to-charge (mz) range(s) for the chromatogram(s). For each
#'     spectrum/retention time, all intensity values within this mz range are
#'     aggregated to result in the intensity value for the spectrum/retention
#'     time. If not specified, the full mz range is considered. See examples
#'     below for details.
#'
#' @param aggregationFun \code{character} defining the function to be used for
#'     intensity value aggregation along the mz dimension. Allowed values are
#'     \code{"sum"} (TIC), \code{"max"} (BPC), \code{"min"} and \code{"mean"}.
#'
#' @param missing \code{numeric(1)} allowing to specify the intensity value for
#'     if for a given retention time (spectrum) no signal was measured within
#'     the mz range. Defaults to \code{NA_real_}.
#'
#' @param BPPARAM Parallelisation backend to be used, which will
#'     depend on the architecture. Default is
#'     \code{BiocParallel::bparam()}.
#'
#' @return \code{chromatogram} returns a \code{\link{Chromatograms}} object with
#'     the number of columns corresponding to the number of files in
#'     \code{object} and number of rows the number of specified ranges (i.e.
#'     number of rows of matrices provided with arguments \code{mz} and/or
#'     \code{rt}).
#' 
#' @author Johannes Rainer
#'
#' @seealso \code{\link{Chromatogram}} and \code{\link{Chromatograms}} for the
#'     classes that represent single and multiple chromatograms.
#'
#' @examples
#' ## Read a test data file.
#' library(msdata)
#' f <- c(system.file("microtofq/MM14.mzML", package = "msdata"),
#'      system.file("microtofq/MM8.mzML", package = "msdata"))
#'
#' ## Read the data as an MSnExp
#' msd <- readMSData(f, msLevel = 1)
#'
#' ## Extract the total ion chromatogram for each file:
#' tic <- chromatogram(msd)
#'
#' tic
#'
#' ## Extract the TIC for the second file:
#' tic[1, 2]
#'
#' ## Plot the TIC for the first file
#' plot(rtime(tic[1, 1]), intensity(tic[1, 1]), type = "l",
#'     xlab = "rtime", ylab = "intensity", main = "TIC")
#'
#' ## Extract chromatograms for a MS data slices defined by retention time
#' ## and mz ranges.
#' rtr <- rbind(c(10, 60), c(280, 300))
#' mzr <- rbind(c(140, 160), c(300, 320))
#' chrs <- chromatogram(msd, rt = rtr, mz = mzr)
#'
#' ## Each row of the returned Chromatograms object corresponds to one mz-rt
#' ## range. The Chromatogram for the first range in the first file is empty,
#' ## because the retention time range is outside of the file's rt range:
#' chrs[1, 1]
#'
#' ## Also the Chromatogram for the second range in the second file is empty
#' chrs[2, 2]
#'
#' ## Get the extracted chromatogram for the first range in the second file
#' chr <- chrs[1, 2]
#' chr
#'
#' plot(rtime(chr), intensity(chr), xlab = "rtime", ylab = "intensity")
setMethod("chromatogram", "MSnExp", function(object, rt, mz,
                                             aggregationFun = "sum",
                                             missing = NA_real_,
                                             BPPARAM = bpparam()){
    if (!missing(rt))
        if (is.null(ncol(rt)))
            rt <- matrix(range(rt), ncol = 2, byrow = TRUE)
    if (!missing(mz))
        if (is.null(ncol(mz)))
            mz <- matrix(range(mz), ncol = 2, byrow = TRUE)
    res <- .extractMultipleChromatograms(object, rt = rt, mz = mz,
                                         aggregationFun = aggregationFun,
                                         missingValue = missing,
                                         BPPARAM = BPPARAM)
    res <- as(res, "Chromatograms")
    res@phenoData <- object@phenoData
    colnames(res) <- rownames(pData(object))
    if (validObject(res))
        res
})
