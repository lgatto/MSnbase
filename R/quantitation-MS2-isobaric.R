quantify_MSnExp <- function(object, method,
                            reporters, strict,
                            BPPARAM,
                            qual,
                            verbose) { ## ignored
    if (any(centroided(object)) & method == "trapezoidation")
        warning("You are quantifying using 'trapezoidation' on centroided data!",
                immediate. = TRUE)
    spectraList <- spectra(object)
    ## Quantification -- creating exprs for assayData slot
    peakData <- bplapply(spectraList, quantify, method, reporters, strict,
                         BPPARAM = BPPARAM)
    .exprs <- do.call(rbind, lapply(peakData, "[[", "peakQuant"))
    rownames(.exprs) <- sub(".peakQuant", "", rownames(.exprs))    
    ## Time consuming - consider removing or caching
    if (qual) {
        qlist <- lapply(peakData, "[[", "curveStats")
        .qual <- do.call(rbind, qlist)
        rownames(.qual) <- sub(".curveStats", "", rownames(.qual))
    } else {
        .qual <- data.frame()
    }
    ## Updating MSnprocess slot
    object@processingData@processing <- c(object@processingData@processing,
                                          paste(reporters@name,
                                                ifelse(strict, " (strict) ", " "),
                                                "quantification by ", method,
                                                ": ", date(), sep=""))
    ## Creating new featureData slot or creating one
    ## if (verbose)
    ##     message("Preparing meta-data")
    fd <- header(object) ## Time consuming - consider caching
    if (nrow(fData(object)) > 0) {
        if (nrow(fData(object)) == length(object)) {
            fd <- combine(fData(object), fd)
        } else {
            warning("Unexpected number of features in featureData slot. Dropping it.")
        }
    }
    ## featureData rows must be reordered to match assayData rows
    .featureData <- new("AnnotatedDataFrame",
                        data = fd[rownames(.exprs), ])
    ## Creating new phenoData slot or creating one
    .phenoData <- new("AnnotatedDataFrame",
                      data = data.frame(mz = reporters@mz,
                          reporters = reporters@name,
                          row.names = reporters@reporterNames))
    if (nrow(pData(object)) > 0) {
        if (nrow(pData(object)) == length(reporters)) {
            .phenoData <- combine(phenoData(object), .phenoData)
        } else {
            ## Here, something more clever should be done, like replicating
            ## old phenoData variables length(reporters) times
            if (verbose)
                message("Original MSnExp and new MSnSet have different number of samples in phenoData. Dropping original.")
        }
    }
    if (verbose)
        message("Creating 'MSnSet' object")
    msnset <- new("MSnSet",
                  qual = data.frame(.qual),
                  exprs = .exprs,
                  experimentData = experimentData(object),
                  phenoData = .phenoData,
                  featureData = .featureData,
                  annotation = "No annotation")

    ## copying processingData
    msnset@processingData <- object@processingData

    ## Updating protocol slot
    if (nrow(protocolData(object)) > 0) {
        if (nrow(protocolData(object)) == length(reporters)) {
            .protocolData <- protocolData(object)
        } else {
            warning("protocolData does not match with reporters. Dropping it.")
        }
    }
    ## Returning shiny MSnSet object
    if (validObject(msnset))
        return(msnset)
}

