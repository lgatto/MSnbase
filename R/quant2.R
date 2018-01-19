setClassUnion("ReporterIonsOrNull", c("ReporterIons", "NULL"))

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

setMethod("show", "QuantitationParam",
          function(object) {
              cat("Quantitation method: '", object@name, "'\n", sep = "")
              cat(" Method:", object@method, "\n")
              if (length(object@msLevel))
                  cat(" MS level:", object@msLevel, "\n")
              if (!is.null(object@reporters))
                  cat(" Reporters:", names(object@reporters), "\n")
          })

IsobaricQuantitation <- function(reporters,
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

SpectralCountingQuantitation <- function(method =
                                             c("count", "SI",
                                               "SIgi", "SIn",
                                               "SAF")) {
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
