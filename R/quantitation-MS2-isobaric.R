quantify_MSnExp <- function(object, method,
                            reporters, strict,
                            BPPARAM,
                            verbose) { ## ignored
    if (any(centroided(object)) & method == "trapezoidation")
        warning("You are quantifying using 'trapezoidation' on centroided data!",
                immediate. = TRUE)
    spectraList <- spectra(object)
    ## Quantification -- creating exprs for assayData slot
    peakData <- bplapply(spectraList, quantify, method, reporters, strict,
                         BPPARAM = BPPARAM)
    .exprs <- do.call(rbind, sapply(peakData, "[", "peakQuant"))
    ## Time consuming - consider removing or caching
    ## .qual <- do.call(rbind, sapply(peakData, "[", "curveStats")) 
    .qual <- data.frame()
    rownames(.exprs) <- sub(".peakQuant", "", rownames(.exprs))
    ## rownames(.qual) <- sub(".curveStats", "", rownames(.qual))
    ## Updating MSnprocess slot
    object@processingData@processing <- c(object@processingData@processing,
                                          paste(reporters@name,
                                                ifelse(strict, " (strict) ", " "),
                                                "quantification by ", method,
                                                ": ", date(), sep=""))
    ## Creating new featureData slot or creating one
    if (verbose)
        message("Preparing meta-data")
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
                  qual = .qual,
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

