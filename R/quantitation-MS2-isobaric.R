quantify_MSnExp <- function(object, method,
                            reporters, strict,
                            parallel, verbose) {
    ## Display progress bar with eapply
    ## TODO - test if using eapply is more efficient in terms of mem/cpu usage
    ## if (verbose) {
    ##   ._cnt <- 1
    ##   pb <- txtProgressBar(min = 0, max = length(object), style = 3)
    ##   ## Quantification -- creating exprs for assayData slot
    ##   peakData <- eapply(assayData(object),function(x) {
    ##     setTxtProgressBar(pb, ._cnt)
    ##     ._cnt <<- ._cnt+1
    ##     quantify(x,method,reporters)
    ##   })
    ##   close(pb)
    ##   rm(pb)
    ##   rm(._cnt)
    ## } else {
    ##   peakData <- eapply(assayData(object),quantify,method,reporters)
    ## }
    detectCores <- registerDoMC <- NULL ## no visible global function definition
    ifelse(verbose, progress <- "text", progress <- "none")
    if (.Platform$OS.type == "windows") {
        parallel <- FALSE
        if (verbose)
            message("Parallel processing not yet supported on Windows.")
    }
    if (any(centroided(object)) & method == "trapezoidation")
        warning("You are quantifying using 'trapezoidation' on centroided data!",
                immediate. = TRUE)
    spectraList <- spectra(object)
    ## Quantification -- creating exprs for assayData slot
    if (length(spectraList) == 1) {
        peakData <- quantify(spectraList[[1]], method, reporters, strict)
        .exprs <- t(peakData$peakQuant)
        rownames(.exprs) <- featureNames(object)
        ## .qual <- as.data.frame(peakData$curveStats)
        ## rownames(.qual) <-
        ##   paste(reporterNames(reporters),
        ##         rep(featureNames(object), each = length(reporters)),
        ##         sep = ".")
        .qual <- data.frame()
    } else {
        if (parallel && require("foreach") && require("doMC") && require("parallel")) {
            registerDoMC(cores = detectCores())
        }
        peakData <- llply(spectraList, quantify, method, reporters, strict,
                          .progress = progress, .parallel = parallel)
        .exprs <- do.call(rbind, sapply(peakData, "[", "peakQuant"))
        ## .qual <- do.call(rbind, sapply(peakData, "[", "curveStats")) ## Time consuming - consider removing or caching
        .qual <- data.frame()
        rownames(.exprs) <- sub(".peakQuant", "", rownames(.exprs))
        ## rownames(.qual) <- sub(".curveStats", "", rownames(.qual))
    }
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
    .featureData <- new("AnnotatedDataFrame", data=fd[rownames(.exprs), ])
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

