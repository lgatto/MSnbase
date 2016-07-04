############################################################
## subsetFeatureDataBy
##
## Convenience function to subset a OnDiskMSnExp featureData based on
## provided subsetting criteria:
## o index: subset by numeric index, logical or character. Character
##    indices are "forwarded" to argument "name". index refers to the
##    rows in the feature data.
## o scanIdx: a numeric is expected, specifying the scan index. If a
##    character vector is provided, it is "forwarded" to
##    argument "name".
## o scanIdxCol: the column of the featureData containing the scan indices.
## o name: a character vector, matching is performed using the row names
##    of the fd data.frame.
## o rtlim: a numeric of length 2 specifying the retention time window
##    from which spectra should be extracted.
## o msLevel: a numeric of MS level.
## Open questions:
## 1) Should we re-calculate the spIdx or the seqNum if we're subsetting?
## Note and comments:
## o At present filtering by multiple criteria (e.g. index and name) is not
##    possible.
subsetFeatureDataBy <- function(fd, index = NULL, scanIdx = NULL,
                                scanIdxCol = "acquisitionNum", name = NULL,
                                rtlim = NULL, msLevel = NULL) {
    valMsg <- validateFeatureDataForOnDiskMSnExp(fd)
    if (is.character(valMsg))
        stop(valMsg)
    ## o index
    if (length(index) > 0) {
        if (is.logical(index)) {
            if (length(index) != nrow(fd))
                stop("If 'index' is a logical vector its length has to match the number of",
                     " rows of the featureData!")
            index <- which(index)
        }
        if (is.numeric(index)) {
            gotIt <- index %in% 1:nrow(fd)
            if (!any(gotIt))
                stop("Provided indices are outside of the allowed range.")
            if (any(!gotIt))
                warning("Some of the provided indices are outside of the allowed range.")
            index <- index[gotIt]
            return(fd[index, , drop = FALSE])
        }
        if (is.character(index))
            name <- index
    }
    ## o scanIdx
    if (length(scanIdx) > 0) {
        if(is.numeric(scanIdx)){
            gotIt <- scanIdx %in% fd[, scanIdxCol]
            if(!any(gotIt))
                stop("None of the provided scan indices are available!")
            if(!all(gotIt))
                warning("Some of the provided scan indices are not available.")
            return(fd[which(fd[, scanIdxCol] %in% scanIdx), , drop=FALSE])
        }
        if (is.character(scanIdx))
            name <- scanIdx
    }
    ## o name: subset by name, match to rownames.
    if (length(name) > 0) {
        gotIt <- name %in% rownames(fd)
        if (!any(gotIt))
            stop("None of the provided names found.")
        if (!all(gotIt))
            warning("Some of the provided names do not match featureData rownames.")
        name <- name[gotIt]
        return(fd[name, , drop=FALSE])
    }
    ## o rtlim: subset by retention time range.
    if (length(rtlim > 0)){
        if (length(rtlim) > 2 | !is.numeric(rtlim))
            stop("Argument 'rtlim' has to be a numeric vector of length 2 specifying",
                 " the retention time window (range).")
        gotIt <- which(fd$retentionTime >= rtlim[1] & fd$retentionTime <= rtlim[2])
        if (length(gotIt) == 0)
            stop("No spectrum within the specified retention time window.")
        fd <- fd[gotIt, , drop=FALSE]
    }
    ## o msLevel
    fd <- fd[fd$msLevel %in% msLevel, , drop = FALSE]
    if (nrow(fd) == 0)
        warning("No spectra with the specified MS level present.")
    return(fd)
}

## Returns either NULL or a character string.
validateFeatureDataForOnDiskMSnExp <- function(x) {
    ## Testing if we've got all the required columns! See issue 105
    ## for a discussion about originalTotIonCurrent and
    ## originalPeaksCount.
    reqCols <- c("fileIdx", "spIdx", "acquisitionNum",
                 "retentionTime", "polarity", "msLevel",
                 "totIonCurrent", "originalPeaksCount",
                 "centroided")
    NotPresent <- reqCols[!(reqCols %in% colnames(x))]
    if (length(NotPresent) > 0)
        return(paste0("Required columns: ",
                      paste(sQuote(NotPresent), collapse = ","),
                      " not present in featureData!"))
    return(NULL)
}
