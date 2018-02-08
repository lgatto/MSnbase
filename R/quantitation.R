## NOTES
##
## Need to focus on quantify,[OnDisk]MSnExp,QuantitationParam.
##
## Quantitation is (currently) MS2 labelled or unlabelled.
##
## Inputs can be a single MS2 spectrum (for isobaric tagging only), or
## a full experiment as MSnExp or OnDiskMSnExp.
##
## Also want to systematically return quantitative features for all
## spectra in the raw data, even when not relevant. There could either
## be an argument in quantify to subset, or this subsetting can be
## done beforehand with filterMsLevel, or afterwards.
## 
## Relevant issues:
## https://github.com/lgatto/MSnbase/issues/130
## https://github.com/lgatto/MSnbase/issues/297
##
## 


##' @rdname quantify-methods
##' @name quantify-methods
##' @aliases quantify
##' 
##' @title Quantifies raw MS data objects
##'
##' @description These methods quantify individual spectra or full MS
##'     experiments. Currently, MS2-level isobar tagging using iTRAQ
##'     and TMT (or any arbitrary peaks of interest, see
##'     [ReporterIons()] and MS2-level label-free quantitation
##'     (spectral counting, spectral index or spectral abundance
##'     factor) are available, as described in
##'     [QuantitationParam()].
##'
##' @details [ReporterIons()] define specific MZ at which peaks are
##'     expected and a window around that MZ value. A peak of interest
##'     is searched for in that window.  Once the range of the curve
##'     is found, quantification is performed.  If no data points are
##'     found in the expected region, `NA` is returned for the
##'     reporter peak MZ.
##'
##'     For spectral counting, spectra that have not been identified
##'     (the corresponding fields in the feature data are populated
##'     with `NA` values) or that have not been uniquely assigned to a
##'     protein (the `nprot` feature variable is greater that 1) are
##'     removed prior to quantitation. The latter does not apply for
##'     method `"count"` (see [QuantitationParam()]) but can be
##'     applied manually with [removeMultipleAssignment()].
##'
##' @param object Raw mass spectrometry data to be quantified. Can be
##'     a single [Spectrum()] or a complete experiment provided as
##'     [MSnExp()] or [OnDiskMSnExp()].
##'
##' @param method An instance of `QuantitationParam` defining the
##'     quantitation method. For isobaric tagging, the method
##'     quantifies peaks defined in [IsobaricTagging()] object. For
##'     MS2-based label-free quantitation, use a [SpectralCounting()]
##'     quantitation parameter.
##'
##'     See [QuantitationParam()] for additional details on how to
##'     define quantiation parameters.
##'
##' @param BPPARAM Support for parallel processing using the
##'     `BiocParallel` infrastructure. When missing (default), the
##'     default registered `BiocParallelParam` parameters are applied
##'     using `bpparam()`. Alternatively, one can pass a valid
##'     `BiocParallelParam` parameter instance: `SnowParam`,
##'     `MulticoreParam`, `DoparParam`, ... see the
##'     `BiocParallel` package for details.
##'
##' @param verbose `logical(1)` defining whether a progress bar should
##'     be displayed when quantifying full experiments. Default is
##'     [isMSnbaseVerbose()].
##'
##' @param ... Further arguments passed to the quantitation functions.
##'
##' @author Laurent Gatto and Sebastian Gibb
##'
##' @return An object of class [MSnSet()] is returned containing the
##'     quantified feature expression and all meta data inherited from
##'     the mass spectrometry experiment data. When quantifying
##'     individual spectra with isobaric tagging, a list of length 2
##'     will be returned. The first element, named `peakQuant`, is a
##'     `numeric` of length equal to the number of reportesr with
##'     quantitation of the reporter peaks.
##'
##'     The second element, named `curveStats`, is a `data.frame`
##'     providing, for each reporter, its curve parameters: maximum
##'     intensity *maxInt*, number of maxima *nMaxInt*, the number of
##'     data points defined the curve *baseLength*, lower and upper MZ
##'     values for the curve (*lowerMz* and *upperMz*), and reporter
##'     reporter and precursor MZ values (when available).
##'
##' @seealso The legacy interface (i.e. that doesn't use
##'     `QuantitationParams`) to `quantify`, in [quantify-legacy].
##' 
##' @md
##' @examples
##' ## On-disk Raw data object
##' library("msdata")
##' f <- proteomics(full.names = TRUE, pattern = "MS3TMT11.mzML")
##' ms <- readMSData(f, mode = "onDisk")
##'
##' ## Isobaric quantitation
##' tmt11 <- IsobaricTagging(reporters = TMT11)
##' msnset <- quantify(ms, tmt11)
##' msnset
##'
##' ## All spectra are present in the quantified data, but only MS3
##' ## spectra have been quantified.
##' head(exprs(msnset))
##' head(exprs(msnset)[msLevel(ms) == 3L, ])
##'
##' ## specifying a custom parallel framework
##' ## bp <- MulticoreParam(2L) # on Linux/OSX
##' ## bp <- SnowParam(2L)      # on Windows
##' ## quantify(ms, tmt11, BPPARAM = bp)
##'
##' ## Quantifying a single spectrum (using the legacy interface only)
##' qty <- quantify(ms[[4]], method = tmt11@method, reporters = tmt11@reporters)
##' qty$curveStats
##' qty$peakQuant
##' ## same as
##' exprs(msnset)[4, ]
##'
##' ## Label-free quantitation
##' ## We need to first add identification data to the raw object
##' data(fdms3tmt11) ## from msdata
##' fData(ms) <- fdms3tmt11
##'
##' ## Spectral counting, setting the peptide sequence feature name to
##' ## 'Sequence' (default is 'sequence'). 
##' sc <- SpectralCounting(pepseq = "Sequence")
##' ## Quantitation for all spectra is returned
##' spc <- quantify(ms, sc)
##' 
##' ## MS2 with a sequence are counted as 1
##' ## MS2 without a sequence are counted as 0
##' ## All other spectra as NA
##' head(exprs(spc))
##' head(fData(ms)[, c("msLevel", "Sequence")])
##'
##' ## Generate a peptide couting
##' ## For example, this peptide was counted 3 times, in spectra 727,
##' ## 733 and 737
##' hvl <- which(fData(spc)$Sequence == "HVLHVQLNRPNK")
##' exprs(spc)[hvl, , drop = FALSE]
##' combineFeatures(spc, fcol = "Sequence", cv = FALSE)
##'
##' ## Total ion count, returns values for all spectra
##' tc <- quantify(ms, SpectralCounting(method = "tic"))
##' head(exprs(tc))
##' 
##' ## It is of course possible to focus only on MS levels of interest
##' ## by either filtering the raw data object or the resulting MSnSet
##' ## using the filterMsLevel method:
##' ms2 <- filterMsLevel(ms, 2L)
##' head(msLevel(ms2))
##' tc2 <- quantify(ms2, SpectralCounting(method = "tic"))
##' head(exprs(tc2))
##' head(exprs(filterMsLevel(tc, 2L)))
##'  
##' ## TODO following ones
##' ## TODO filter MS2 and document this
##' # si <- quantify(msexp, method = "SIn")
##' # processingData(si)
##' # exprs(si)
##'
##' # saf <- quantify(msexp, method = "NSAF")
##' # processingData(saf)
##' # exprs(saf)
NULL
##> NULL


## What about the former qual argument: Should the \code{qual} slot be
## populated. Default is \code{TRUE}.

quantify2 <- function(object,
                      params,
                      BPPARAM,
                      verbose,
                      ...) {
    stopifnot(inherits(params, "QuantitationParam"))
    if (missing(BPPARAM)) {
        BPPARAM <- bpparam()
        if (verbose) message("Using default parallel backend: ",
                             class(BPPARAM)[1])
    }
    if (params@name == "IsobaricTagging") {
        if (!length(params@msLevel))
            params@msLevel <- max(msLevel(object))
        obj2 <- filterMsLevel(object, params@msLevel)        
        if (params@method == "max") {
            if (!verbose)
                suppressMessages(e <- quantify_OnDiskMSnExp_max(obj2,
                                                                params@reporters,
                                                                params@wd,
                                                                BPPARAM))
            else e <- quantify_OnDiskMSnExp_max(obj2, params@reporters,
                                                params@wd, BPPARAM)
        } else {
            quantify_MSnExp(obj2,
                            params@method,
                            params@reporters,
                            params@strict,
                            BPPARAM, qual = FALSE,
                            verbose)
        }
        ans <- matrix(NA_real_,
                      nrow = length(object),
                      ncol = ncol(e),
                      dimnames = list(featureNames(object),
                                      sampleNames(e))) 
        ans[featureNames(e), ] <- exprs(e)
        ans <- MSnSet(exprs = ans,
                      fData = fData(object),
                      pData = pData(e))
        ans@processingData <- e@processingData    
        return(ans)
    } else if (params@name == "SpectralCounting") {
        if (params@method == "count") {
            count_MSnSet(object,
                         pepseq = params@pepseq,
                         removeNoId. = FALSE)
        } else if (params@method == "tic") {
            tic_MSnSet(object)

        } else {
            ## the following assumes that the appropriate fcols
            ## are available in the QuantitationParam object
            object <- utils.removeNoIdAndMultipleAssignments(object)
            if (params@method %in% c("SI", "SIgi", "SIn"))
                SI(object, params@method,
                   groupBy = params@dbaccess,
                   plength = params@plength,
                   verbose = verbose)
            else SAF(object, params@method,
                     groupBy = params@dbaccess,
                     plength = params@plength,
                     verbose = verbose)
        }
    } else if (params@name == "LFQ") {
        stop("LFQ currently not implemented.")
    } else
        stop("Quantitation method not recognised.")
}
