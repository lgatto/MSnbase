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
          function(object, method, ...) quantify(object, method, ...))

transferQuantToPrecursorScanNum <- function(x) {
    e <- matrix(NA_real_, ncol = ncol(x), nrow = nrow(x))
    rownames(e) <- fData(x)[, "acquisitionNum"]
    ms3 <- fData(x)$msLevel == 3L
    e[as.character(fData(x)$precursorScanNum[ms3]), ] <- exprs(x)[ms3, ]
    e
}
