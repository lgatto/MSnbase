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
##'     [QuantificationParam()].
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
##'     label-free quantitation, defined by a [SpectralCounting()]
##'     quantitation parameter, the respective quantitation methods
##'     and normalisations are applied to the spectra. These methods
##'     require two additional arguments passed in `...`, namely the
##'     protein accession of identifiers (`fcol`, with detault value
##'     `"DatabaseAccess"`) and the protein lengths (`plength`, with
##'     default value `"DBseqLength"`). These values are available of
##'     the identification data had been collated using
##'     [addIdentificationData()]. 
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
##' @param pepseq `character(1)` giving the peptide sequence column in
##'     the feature data. Default is `"sequence"`.
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
##' @seealso The legacy interface to `quantify` is still availabel in
##'     [quantify-legacy].
##' 
##' @md
##' @examples
##' ## Quantifying a full experiment using iTRAQ4-plex tagging
##' data(itraqdata)
##' msnset <- quantify(itraqdata, method = "trap", reporters = iTRAQ4)
##' msnset
##'
##' ## specifying a custom parallel framework
##' ## bp <- MulticoreParam(2L) # on Linux/OSX
##' ## bp <- SnowParam(2L) # on Windows
##' ## quantify(itraqdata[1:10], method = "trap", iTRAQ4, BPPARAM = bp)
##'
##' ## Checking for non-quantified peaks
##' sum(is.na(exprs(msnset)))
##'
##' ## Quantifying a single spectrum
##' qty <- quantify(itraqdata[[1]], method = "trap", iTRAQ4[1])
##' qty$peakQuant
##' qty$curveStats
##'
##' ## Label-free quantitation
##' ## Raw (mzXML) and identification (mzid) files
##' quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
##'                  full.name = TRUE, pattern = "mzXML$")
##' identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
##'                  full.name = TRUE, pattern = "dummyiTRAQ.mzid")
##' msexp <- readMSData(quantFile)
##' msexp <- addIdentificationData(msexp, identFile)
##' fData(msexp)$DatabaseAccess
##'
##' si <- quantify(msexp, method = "SIn")
##' processingData(si)
##' exprs(si)
##'
##' saf <- quantify(msexp, method = "NSAF")
##' processingData(saf)
##' exprs(saf)
NULL
##> NULL


## What about the former qual argument: Should the \code{qual} slot be
## populated. Default is \code{TRUE}.

quantify2 <- function(object,
                      params,
                      BPPARAM,
                      pepseq,
                      verbose,
                      ...) {
    stopifnot(inherits(params, "QuantitationParam"))
    if (missing(BPPARAM)) {
        BPPARAM <- bpparam()
        if (verbose) message("Using default parallel backend: ",
                             class(BPPARAM)[1])
    }
    if (params@name == "IsobaricTagging") {
        if (params@method != "max")
            stop("Not yet implemented - see issue #130")
        if (!length(params@msLevel))
            params@msLevel <- max(msLevel(object))
        obj2 <- filterMsLevel(object, params@msLevel)        
        if (!verbose)
            suppressMessages(e <- quantify_OnDiskMSnExp_max(obj2,
                                                            params@reporters,
                                                            params@wd,
                                                            BPPARAM))
        else e <- quantify_OnDiskMSnExp_max(obj2, params@reporters,
                                            params@wd, BPPARAM)
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
        stop("TODO")
    } else if (params@name == "LFQ") {
        stop("LFQ currently not implemented.")
    } else
        stop("Quantitation method not recognised.")
}
