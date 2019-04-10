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
                  processingData = processingData(object),
                  annotation = "No annotation")

    ## Updating protocol slot
    if (nrow(protocolData(object)) > 0) {
        if (nrow(protocolData(object)) == length(reporters)) {
            msnset@protocolData <- protocolData(object)
        } else {
            warning("protocolData does not match with reporters. Dropping it.")
        }
    }
    ## Returning shiny MSnSet object
    if (validObject(msnset))
        return(msnset)
}


## quantifies highest peaks pk in regions [pks - wd, pks + wd] in
## spectra spi in file f using max of peak (irrespective if spectrum
## is centroided or profile mode)
fastquant_max <- function(f, pk, spi, wd = 0.5) {
    ramp <- .openMSfile(f)
    on.exit(close(ramp))
    pks <- peaks(ramp, spi)
    if (length(spi) == 1)
        pks <- list(pks)
    mzs <- res <- matrix(NA_real_,
                         ncol = length(pk),
                         nrow = length(pks))
    for (i in seq_along(pks)) {
        for (ii in seq_along(pk)) {
            k <- abs(pks[[i]][, 1L] - pk[ii]) <= wd
            if (any(k)) {
                j <- which.max(pks[[i]][k, 2])
                res[i, ii] <- pks[[i]][k, , drop = FALSE][j, 2]
                mzs[i, ii] <- pks[[i]][k, , drop = FALSE][j, 1]
            }
        }
    }
    list(exprs = res,
         mzs = mzs)
}


## fastquant_max_1 <- function(f, pk, spi, wd = 0.5) {
##     ramp <- openMSfile(f)
##     on.exit(close(ramp))
##     pks <- peaks(ramp, spi)
##     if (length(spi) == 1)
##         pks <- list(pks)
##     res <- matrix(NA_real_,
##                   ncol = length(pk),
##                   nrow = length(pks))
##     for (i in seq_along(pks)) {
##         for (ii in seq_along(pk)) {
##             k <- pks[[i]][, 1] >= pk[ii] - wd & pks[[i]][, 1] <= pk[ii] + wd
##             if (any(k)) {
##                 res[i, ii] <- max(pks[[i]][k, 2])
##             }
##         }
##     }
##     res
## }

quantify_OnDiskMSnExp_max <- function(object, reporters,
                                      wd,
                                      BPPARAM) {
    if (missing(reporters) | !inherits(reporters, "ReporterIons"))
        stop("Valid reporters required.")
    if (missing(wd)) wd <- width(reporters)
    fDataPerFile <- split(fData(object), f = fData(object)$fileIdx)
    if (missing(BPPARAM)) BPPARAM <- bpparam()

    res <- bplapply(fDataPerFile,
                    FUN = function(fdf)
                        fastquant_max(
                            f = fileNames(object)[fdf$fileIdx[1]],
                            pk = mz(reporters),
                            spi = fdf$spIdx,
                            wd),
                    BPPARAM = BPPARAM)
    ## quantitation data
    e <- do.call(rbind, lapply(res, "[[", "exprs"))
    colnames(e) <- reporterNames(reporters)
    rownames(e) <- unlist(lapply(fDataPerFile, rownames),
                          use.names = FALSE)
    e <- e[featureNames(object), ]
    ## MZ at max (for reporter accuracy QC)
    mzs <- do.call(rbind, lapply(res, "[[", "mzs"))
    colnames(mzs) <- reporterNames(reporters)
    rownames(mzs) <- unlist(lapply(fDataPerFile, rownames),
                            use.names = FALSE)
    mzs <- mzs[featureNames(object), ]
    rm(res)

    .phenoData <- new("AnnotatedDataFrame",
                      data = data.frame(mz = mz(reporters),
                                        reporters = names(reporters),
                                        row.names = reporterNames(reporters)))

    ## This actually fails if the number of files (rows) in the MSnExp
    ## matches the number of reporter ions in the MSnSet, as it tried
    ## to combine an NAnnotatedDataFrame (now deprecated) and a
    ## AnnotatedDataFrame.
    ##
    ## if (nrow(pData(object)) > 0) {
    ##     if (nrow(pData(object)) == length(reporters)) {
    ##         .phenoData <- combine(phenoData(object), .phenoData)
    ##     } else {
    ##         ## Here, something more clever should be done, like replicating
    ##         ## old phenoData variables length(reporters) times
    ##         msg <- paste(strwrap(paste0("Original MSnExp and new MSnSet have ",
    ##                                     "different number of samples in ",
    ##                                     "phenoData. Dropping original.")),
    ##                      collapse = "\n")
    ##         message(msg)
    ##     }
    ## }

    ans <- new("MSnSet",
               exprs = e,
               featureData = featureData(object),
               phenoData = .phenoData,
               processingData = processingData(object))
    fData(ans)$reporterMzs <- mzs
    ans <- logging(ans, paste0("Fast ", names(reporters),
                               " quantitation by max"))
    if (validObject(ans)) ans
}
