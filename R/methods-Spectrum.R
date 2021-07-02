##################################################################
## Methods for Spectrum class and children
setMethod("initialize",
          "Spectrum",
          function(.Object, ..., mz, intensity, peaksCount) {
              ## Add Spectrum class version - set in Spectrum1/2 to
              ## keep suber/sub class order
              ## classVersion(.Object)["Spectrum"] <- getClassVersion("Spectrum")
              if (xor(!missing(mz), !missing(intensity)))
                  stop("'mz' and 'intensity' or none required.")
              if (!missing(mz) & !missing(intensity)) {
                  ## issue #180: use radix sorting for R >= 3.3
                  ## o <- base::order(mz, method = "radix")
                  ## o <- base::order(mz, method = MSnbaseOptions()$sortMethod)
                  ## To fix it for sure we're dropping the method parameter
                  ## completely
                  o <- base::order(mz)
                  .Object <- callNextMethod(.Object,
                                            ...,
                                            mz = mz[o],
                                            intensity = intensity[o],
                                            peaksCount = length(mz))
              } else .Object <- callNextMethod(.Object, ...)
              if (.Object@tic == 0)
                  .Object@tic <- sum(.Object@intensity)
              if (validObject(.Object))
                  .Object
          })

setMethod("show", "Spectrum",
          function(object) {
              if (msLevel(object) == 1) show_Spectrum1(object)
              else show_Spectrum2(object)
              invisible(NULL)
          })

setMethod("plot", c("Spectrum", "missing"),
          function(x, y, ...) {
              if (msLevel(x) == 1) plot_Spectrum1(x, ...)
              else plot_Spectrum2(x, ...)
          })

setMethod("plot", c("Spectrum", "Spectrum"),
          function(x, y, ...) {
              plotSpectrumVsSpectrum(list(x, y), ...)
          })

setMethod("plot", c("Spectrum2", "character"),
          function(x, y, ...) {
              .plotSingleSpectrum(x, y, ...)
          })

setMethod("clean",
          signature = signature("Spectrum"),
          function(object, all = FALSE, msLevel.)
              clean_Spectrum(object, all, msLevel. = msLevel.))

setMethod("removePeaks", "Spectrum",
          function(object, t, msLevel.)
              removePeaks_Spectrum(object, t, msLevel.))

setMethod("precursorMz", "Spectrum",
          function(object) {
              if (msLevel(object) > 1)
                  return(object@precursorMz)
              stop("No precursor MZ value for MS1 spectra.")
          })

setMethod("precursorCharge", "Spectrum",
          function(object) {
              if (msLevel(object) > 1)
                  return(object@precursorCharge)
              stop("No precursor charge value for MS1 spectra.")
          })

setMethod("precursorIntensity", "Spectrum",
          function(object) {
              if (msLevel(object) > 1)
                  return(object@precursorIntensity)
              stop("No precursor data for MS1 spectra.")
          })

setMethod("acquisitionNum", "Spectrum",
          function(object) object@acquisitionNum)
setMethod("scanIndex", "Spectrum",
          function(object) object@scanIndex)

setMethod("precScanNum", "Spectrum",
          function(object) {
              if (msLevel(object) > 1)
                  return(object@precScanNum)
              stop("This is already an MS1 spectrum.")
          })
setMethod("precAcquisitionNum", "Spectrum",
          function(object) precScanNum(object))

setMethod("rtime", "Spectrum", function(object) object@rt)

setMethod("peaksCount",
          signature("Spectrum", "missing"),
          function(object) object@peaksCount)

setMethod("msLevel", "Spectrum", function(object) object@msLevel)
setMethod("collisionEnergy", "Spectrum",
          function(object) {
              if (msLevel(object) > 1)
                  return(object@collisionEnergy)
              stop("No collision energy for MS1 spectra.")
          })
setMethod("intensity", "Spectrum", function(object) object@intensity)
setMethod("mz", "Spectrum", function(object) object@mz)

setMethod("tic", "Spectrum", function(object) object@tic)

setMethod("ionCount", "Spectrum", function(object) sum(object@intensity))

setMethod("trimMz",
          signature = signature("Spectrum", "numeric"),
          function(object, mzlim, ...)
              trimMz_Spectrum(object, mzlim, ...))

setMethod("filterMz", "Spectrum",
          function(object, mz, msLevel., ...) {
              return(trimMz_Spectrum(object, mzlim = mz,
                                     msLevel. = msLevel., ...))
          })

setMethod("quantify",
          signature = signature("Spectrum"),
          function(object,
                   method = c("trapezoidation", "max", "sum"),
                   reporters,
                   strict = FALSE) {
              if (!inherits(reporters, "ReporterIons"))
                  stop("Argument 'reporters' must inherit from 'ReporterIons' class.")
              quantify_Spectrum(object, match.arg(method), reporters, strict)
          })

## setMethod("curveStats","Spectrum",
##           function(object,reporters) curveStats_Spectrum(object,reporters))

setReplaceMethod("precursorCharge",
                 signature(object = "Spectrum",
                           value = "integer"),
                 function(object, value) {
                     object@precursorCharge <- value
                     if (validObject(object))
                         return(object)
                 })

setMethod("fromFile", "Spectrum", function(object) object@fromFile)

setMethod("polarity", "Spectrum",
          function(object) return(object@polarity))

setAs("Spectrum", "data.frame",
      function (from)
          data.frame(mz = mz(from),
                     i = intensity(from))
      )

as.data.frame.Spectrum <- function(x, row.names=NULL, optional=FALSE, ...)
    as(x, "data.frame")

setMethod("centroided", "Spectrum",
          function(object, na.fail = FALSE) {
              if (na.fail & is.na(object@centroided))
                  stop("Mode is undefined. See ?isCentroided for details.",
                       call. = FALSE)
              object@centroided
          })

setReplaceMethod("centroided",
                 signature(object = "Spectrum",
                           value = "logical"),
                 function(object, value) {
                     object@centroided <- value
                     if (validObject(object))
                         return(object)
                 })

setMethod("smoothed", "Spectrum", function(object) object@smoothed)
setReplaceMethod("smoothed",
                 signature(object = "Spectrum",
                           value = "logical"),
                 function(object, value) {
                     object@smoothed <- value
                     if (validObject(object))
                         return(object)
                 })


setMethod("normalize", "Spectrum",
          function(object, method = c("max", "sum"), ...) {
              normalise_Spectrum(object, method = match.arg(method))
          })

setMethod("normalize", "Spectrum2",
          function(object,
                   method = c("max", "sum", "precursor"),
                   precursorIntensity,
                   ...) {
              method <- match.arg(method)
              if (method == "precursor") {
                  precursorIntensity <- ifelse(missing(precursorIntensity),
                                               object@precursorIntensity,
                                               precursorIntensity)
                  return(normalise_Spectrum(object,
                                            method = "value",
                                            value = precursorIntensity))
              } else {
                  return(callNextMethod(object, method, ...))
              }
          })

normalise <- normalize

setMethod("bin", "Spectrum",
          function(x, binSize = 1L,
                   breaks = seq(floor(min(mz(x))),
                                ceiling(max(mz(x))),
                                by = binSize),
                   msLevel.) {
              bin_Spectrum(x, binSize = binSize, breaks = breaks,
                           msLevel. = msLevel.)
          })

setMethod("calculateFragments", c("character", "Spectrum2"),
          function(sequence, object, ...) {
              calculateFragments_Spectrum2(sequence, object, ...)
          })

setMethod("compareSpectra", c("Spectrum", "Spectrum"),
          function(x, y, fun=c("common", "cor", "dotproduct"),
                   ...) {
              compare_Spectra(x, y, fun = fun, ...)
          })

setMethod("estimateNoise", "Spectrum",
          function(object, method = c("MAD", "SuperSmoother"), ...) {
              estimateNoise_Spectrum(object, method = match.arg(method), ...)
          })

setMethod("pickPeaks", "Spectrum",
          function(object, halfWindowSize = 3L,
                   method = c("MAD", "SuperSmoother"),
                   SNR = 0L, refineMz = c("none", "kNeighbors",
                                          "kNeighbours", "descendPeak"),
                   msLevel. = msLevel(object), ...) {
              pickPeaks_Spectrum(object, halfWindowSize = halfWindowSize,
                                 method = match.arg(method), SNR = SNR,
                                 refineMz = match.arg(refineMz),
                                 msLevel. = msLevel., ...)
          })

setMethod("smooth", "Spectrum",
          function(x, method = c("SavitzkyGolay", "MovingAverage"),
                   halfWindowSize = 2L, msLevel. = msLevel(x),
                   ...) {
              smooth_Spectrum(x, method = match.arg(method),
                              halfWindowSize = halfWindowSize,
                              msLevel. = msLevel., ...)
          })

setMethod("removeReporters", "Spectrum",
          function(object, reporters = NULL, clean = FALSE, ...) {
              if (msLevel(object) > 1)
                  return(removeReporters_Spectrum2(object, reporters, clean))
              ## stop("No reporters to remove for MS1 spectra.")
              ## Instead of stopping we show a (conditional) warning
              ## See also issue #161
              dots <- list(...)
              if (!((length(dots$suppressWarnings) > 0) && dots$suppressWarnings))
                  warning("No reporters to remove for MS1 spectra.")
              ## Return the Spectrum as-is for MS1
              return(object)
          })

setMethod("isEmpty", "Spectrum",
          function(x) length(x@mz) == 0)

setMethod("isCentroided", "Spectrum",
          function(object, ...)
              .isCentroided(as(object, "data.frame"), ...))

#' @title Estimate the m/z resolution of a spectrum
#'
#' @aliases estimateMzResolution
#'
#' @description
#'
#' `estimateMzResolution` estimates the m/z resolution of a profile-mode
#' `Spectrum` (or of all spectra in an [MSnExp] or [OnDiskMSnExp] object.
#' The m/z resolution is defined as the most frequent difference between a
#' spectrum's m/z values.
#'
#' @note
#'
#' This assumes the data to be in profile mode and does not return meaningful
#' results for centroided data.
#'
#' The estimated m/z resolution depends on the number of ions detected in a
#' spectrum, as some instrument don't measure (or report) signal if below a
#' certain threshold.
#'
#' @param object either a `Spectrum`, `MSnExp` or `OnDiskMSnExp` object.
#'
#' @param ... currently not used.
#'
#' @return `numeric(1)` with the m/z resolution. If called on a `MSnExp` or
#' `OnDiskMSnExp` a `list` of m/z resolutions are returned (one for
#' each spectrum).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @rdname estimateMzResolution
#'
#' @examples
#'
#' ## Load a profile mode example file
#' library(BiocParallel)
#' register(SerialParam())
#' library(msdata)
#' f <- proteomics(full.names = TRUE,
#'     pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
#'
#' od <- readMSData(f, mode = "onDisk")
#'
#' ## Estimate the m/z resolution on the 3rd spectrum.
#' estimateMzResolution(od[[3]])
#'
#' ## Estimate the m/z resolution for each spectrum
#' mzr <- estimateMzResolution(od)
#'
#' ## plot the distribution of estimated m/z resolutions. The bimodal
#' ## distribution represents the m/z resolution of the MS1 (first peak) and
#' ## MS2 spectra (second peak).
#' plot(density(unlist(mzr)))
setMethod("estimateMzResolution", "Spectrum", function(object, ...) {
    .estimate_mz_resolution(object@mz)
})
