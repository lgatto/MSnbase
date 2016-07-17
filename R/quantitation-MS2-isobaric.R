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


## quantifies peaks pk in regions [pks - wd, pks + wd] at half window
## size hws in spectra spi in file f using max of peak (irrespective
## if spectrum is centroided or profile mode)
fastquant_max <- function(f, pk, spi, hws, wd = 0.5) {
    ramp <- openMSfile(f)
    on.exit(close(ramp))
    pks <- peaks(ramp, spi)
    res <- matrix(NA_real_,
                  ncol = length(pk),
                  nrow = length(pks))
    for (i in seq_along(pks)) {
        for (ii in seq_along(pk)) {
            k <- pks[[i]][, 1] >= pk[ii] - wd & pks[[i]][, 1] <= pk[ii] + wd
            if (any(k)) {
                .pks <- pks[[i]][k, , drop = FALSE]
                mx <- MALDIquant:::.localMaxima(.pks[, 2], hws)
                res[i, ii] <- .pks[mx, 2]
            }
        }
    }
    res
}


quantify_OnDiskMSnExp_max <- function(object, reporters,
                                      hws = 20L, wd = 0.05,
                                      BPPARAM) {
    if (missing(reporters) | !inherits(reporters, "ReporterIons"))
        stop("Valid reporters required.")
    fDataPerFile <- split(fData(object), f = fData(object)$fileIdx)
    if (missing(BPPARAM)) BPPARAM <- bpparam()

    res <- bplapply(fDataPerFile,
                    FUN = function(fdf)
                        fastquant_max(
                            f = fileNames(object)[fdf$fileIdx[1]],
                            pk = mz(reporters),
                            spi = fdf$spIdx,
                            hws, wd),
                    BPPARAM = BPPARAM)
    res <- do.call(rbind, res)
    colnames(res) <- reporterNames(reporters)
    rownames(res) <- unlist(lapply(fDataPerFile, rownames), use.names = FALSE)
    res <- res[featureNames(object), ]

    .phenoData <- new("AnnotatedDataFrame",
                      data = data.frame(mz = mz(reporters),
                                        reporters = names(reporters),
                                        row.names = reporterNames(reporters)))

    if (nrow(pData(object)) > 0) {
        if (nrow(pData(object)) == length(reporters)) {
            .phenoData <- combine(phenoData(object), .phenoData)
        } else {
            ## Here, something more clever should be done, like replicating
            ## old phenoData variables length(reporters) times
            message(paste(strwrap(paste0("Original MSnExp and new MSnSet have ",
                                         "different number of samples in phenoData. ",
                                         "Dropping original.")), collapse = "\n"))
        }
    }

    ans <- new("MSnSet", exprs = res,
               featureData = featureData(object),
               phenoData = .phenoData)
    ans <- MSnbase:::logging(ans, paste0("Fast ", names(reporters),
                                         " quantitation by max"))
    if (validObject(ans)) ans
}
