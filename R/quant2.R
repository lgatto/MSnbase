setClassUnion("ReporterIonsOrNull", c("ReporterIons", "NULL"))

##' @aliases QuantitationParam IsobaricTagging SpectralCounting
##' 
##' @title Quantitation parameter class
##' 
##' @description `QuantitationParam` object are used to define the
##'     type of quantitation to perform. These objects are then passed
##'     with the raw data object to the [quantify()]
##'     method. `QuantitationParam` objects are rarely created
##'     directly; instead, users should use the predefined constructor
##'     functions `IsobaricTagging` and `SpectralCounting`.
##'
##' @seealso quantify
##' 
##' @slot msLevel `integer(1)` setting the MS level at which
##'     quantitation should be done. It is generally set automatically
##'     for a method (for instance 2 for spectral counting) or from
##'     the data, if absent in the `QuantitationParam` object (the
##'     highest MS level in the data, i.e. 2 or 3, for isobaric
##'     tagging).
##' 
##' @slot reporters The reporter ions to use, as an
##'     [ReporterIons()]. For isobaric tagging quantation only.
##' 
##' @slot method `character(1)` defining the method of
##'     quantitation. For isobaric tagging, this would be one of
##'     `"max"` (default), `"trapezoidation"`, or `"sum"`. For
##'     spectral counting, one of `"count"` (default), `"SI"`,
##'     `"SIgi"`, `"SIn"`, `"SAF"` or `"NSAF"`.
##' 
##' @slot methargs A `list` of additional argument applied to the
##'     quantitation method.
##' 
##' @slot name `characater(1)` naming the quantiation parameter.
##' 
##' @slot wd `numeric(1)` defining the width around the experter
##'     isobaric tag mass to seach for a peak. If omitted, is
##'     extracted from the `ReporterIons` object. For isobaric tagging
##'     quantation only.
##' 
##' @slot strict `logical(1)` defining if a peak should be quantified
##'     beyond the reporter tag +/- the width. For isobaric tagging
##'     quantation only.
##' 
##' @slot .__classVersion__ The version of the `QuantitationParam`
##'     class definition.
##' 
##' @author Laurent Gatto
##' 
##' @rdname QuantitationParam-class
##' 
##' @md
##' 
##' @examples
##' ## Isobaric tagging using iTRAQ 4-plex
##' IsobaricTagging(iTRAQ4)
##' ## Isobaric tagging using TMT 11-plex at the MS3 level
##' IsobaricTagging(TMT11, msLevel = 3L)
##' ## Isobaric tagging using TMT 11-plex at the MS2 level
##' IsobaricTagging(TMT11, msLevel = 2L)
##' 
##' ## Spectral counting, raw counts
##' SpectralCounting()
##' ## Spectral counting, normalised spectral abundance factor
##' SpectralCounting("NSAF")
.QuantitationParam <- 
    setClass("QuantitationParam",
             slots = c(msLevel = "integer",
                       reporters = "ReporterIonsOrNull",
                       method = "character",
                       methargs = "list",
                       name = "character",
                       wd = "numeric",
                       strict = "logical"),
             contains = "Versioned",
             prototype = prototype(
                 name = "Quantitation method",
                 new("Versioned",
                     versions = c(classVersion("ReporterIons"),
                                  QuantitationParam = "0.1.0"))
             ))

##' @rdname QuantitationParam-class
setMethod("show", "QuantitationParam",
          function(object) {
              cat("Quantitation method: '", object@name, "'\n", sep = "")
              cat(" Method:", object@method, "\n")
              if (length(object@msLevel))
                  cat(" MS level:", object@msLevel, "\n")
              if (!is.null(object@reporters))
                  cat(" Reporters:", names(object@reporters), "\n")
          })

##' @rdname QuantitationParam-class
IsobaricTagging <- function(reporters,
                            msLevel, 
                            method =
                                c("max",
                                  "trapezoidation",
                                  "sum")) {
    if (missing(reporters))
        stop("Please provide the isobaric reporters.")
    if (missing(msLevel)) {
        if (isMSnbaseVerbose())
            message("No MS level provided. Will be determined from data.")
        msLevel <- integer()
    }
    method <- match.arg(method)
    .QuantitationParam(msLevel = msLevel,
                       reporters = reporters,
                       method = method,
                       wd = width(reporters),
                       strict = FALSE,
                       name = "IsobaricTagging")
}

##' @rdname QuantitationParam-class
SpectralCounting <- function(method =
                                 c("count", "SI",
                                   "SIgi", "SIn",
                                   "SAF", "NSAF")) {
    method <- match.arg(method)
    .QuantitationParam(msLevel = 2L,
                       method = method,
                       name = "SpectralCounting")
}

quantify2 <- function(object,
                      params,
                      BPPARAM,
                      verbose = isMSnbaseVerbose(),
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


setMethod("quantify", c("OnDiskMSnExp", "QuantitationParam"),
          function(object, method, ...) quanitfy2(object, method, ...))

transferQuantToPrecursorScanNum <- function(x) {
    e <- matrix(NA_real_, ncol = ncol(x), nrow = nrow(x))
    rownames(e) <- fData(x)[, "acquisitionNum"]
    ms3 <- fData(x)$msLevel == 3L
    e[as.character(fData(x)$precursorScanNum[ms3]), ] <- exprs(x)[ms3, ]
    e
}
