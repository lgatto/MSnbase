## Read chromatogram data from an mzML.
## Problem is that the function to read the data is only implemented for pwiz.

## Read SRM data
## Requires mzML files with chromatogram data recorded in selected reaction
## monitoring mode. TIC chromatograms are not loaded.
## SRM: Q1, Q3, start end (acquisition time).
readSRMData <- function(files, pdata = NULL,
                        verbose = isMSnbaseVerbose()) {
    files <- normalizePath(files)
    ## Read first the header of all files.
    hdr_list <- lapply(files, function(x) {
        msf <- MSnbase:::.openMSfile(x)
        if (!is(msf, "mzRpwiz"))
            stop("Can only extract chromatogram information from a mzML file",
                 " using the 'proteowizard' backend")
        hdr <- chromatogramHeader(msf)
        close(msf)
        hdr[!is.na(hdr$precursorIsolationWindowTargetMZ) |
            !is.na(hdr$productIsolationWindowTargetMZ), , drop = FALSE]
    })
    ## Check if we've got any file without SRM data:
    lens <- unlist(lapply(hdr_list, nrow))
    if (any(lens == 0))
        stop("file(s) ", paste0("'", files[lens == 0], "'", collapse = ", "),
             " do not contain SRM chromatogram data")
    ## Now combine the individual headers into a feature data data.frame
    fdata <- .combine_data.frame(hdr_list,
                                 cols = c("precursorIsolationWindowTargetMZ",
                                          "productIsolationWindowTargetMZ"))


    
    ## For SRM data we want to group chromatograms by polarity, precursorMZ,
    ## productMZ
    f_data <- unique(do.call(
        rbind, hdr_list)[, c("polarity",
                             "precursorIsolationWindowTargetMZ",
                             "productIsolationWindowTargetMZ")])
    f_ids <- paste0(.polarity_char(f_data$polarity),
                    " Q1=", f_data$precursorIsolationWindowTargetMZ,
                    " Q3=", f_data$productIsolationWindowTargetMZ)
    ## Loop again through the files to read the data and return Chromatograms
    ## in the order matching the f_data.
    ## Need to get a list() of Chromatogram objects. These are then passed to
    ## the Chromatograms function along with ncol.
    chrs <- mapply(
        files, hdr_list, 1:length(files),
        FUN = function(x, hdr, idx) {
            ## It's not SRM data if precursorIsolationWindowTargetMZ and
            ## productIsolationWindowTargetMS is both NA!
            ## Get the index of the chromatograms in the feature data.
            current_f_ids <- paste0(.polarity_char(hdr$polarity),
                                    " Q1=",
                                    hdr$precursorIsolationWindowTargetMZ,
                                    " Q3=",
                                    hdr$productIsolationWindowTargetMZ)
            ## If these are not unique we've got a problem, i.e. the data is
            ## most likely not SRM data.
            if (length(current_f_ids) != length(unique(current_f_ids)))
                stop("Chromatogram data in file ", basename(x), " does not ",
                     "appear to represent SRM data")

            ## AAAAAH, got Q1=180.1 Q3=162.996 twice! with slightly different
            ## retention time windows! Would have to add two rows for this guy!
            ## Load the data.
            msf <- MSnbase:::.openMSfile(x)
            chr <- mapply(
                chromatogram(msf), split.data.frame(hdr, 1:nrow(hdr)),
                FUN = function(dat, hd) {
                    Chromatogram(
                        rtime = dat[, 1],
                        intensity = dat[, 2],
                        precursorMz = hd$precursorIsolationWindowTargetMZ,
                        productMz = hd$precursorIsolationWindowTargetMZ,
                        fromFile = idx)
                })
            close(msf)
        })
}

#' @title Combine `data.frame`s without collapsing non-unique rows
#'
#' @description
#'
#' Combine a list of `data.frame`s into a single `data.frame`. In
#' contrast to a *unique* combination in which each unique row will only be
#' present once, this function does not reduce redundant rows to a single row,
#' but returns them. In other words, if three rows in one of the input
#' `data.frame` contain the same values, the result `data.frame` also has the
#' same 3 rows, while `unique` would collapse them into a single row. If rows
#' with the same values are present multiple times in multiple input
#' `data.frame`'s, e.g. 3 times in one `data.frame` and 2 times in
#' another one, the result `data.frame` will contain the rows the maximum time
#' they are present across tables. See examples for details.
#'
#' @note
#'
#' The ordering of a combined `data.frame` from `data.frame`s with some or
#' all columns being `factors` might not be as expected. It is thus suggested
#' to avoid passing `data.frame`s with `factors` to the function.
#' 
#' @param x
#'
#' @param cols `character` defining the columns that should be returned.
#'
#' @return `data.frame`
#' @noRd
#'
#' @md
#' 
#' @author Johannes Rainer
#'
#' @examples
#'
#' A <- data.frame(a = c(1, 1, 2, 3, 4, 5, 6), b = c(2, 2, 3, 4, 5, 5, 6))
#' B <- data.frame(a = c(1, 3, 3, 3, 3, 4, 5), b = c(2, 4, 4, 4, 5, 5, 5))
#' C <- data.frame(a = c(3, 3, 4, 4), b = c(4, 5, 5, 5))
#'
#' x <- list(A, B, C)
#'
#' .combine_data.frame(x)
#'
#' D <- data.frame(a = c("z", "b", "a", "g", "g"),
#'                 b = c(1, 2, 3, 4, 4), stringsAsFactors = FALSE)
#' E <- data.frame(a = c("g", "a", "d", "g"),
#'                 b = c(4, 3, 1, 1), stringsAsFactors = FALSE)
#' .combine_data.frame(list(D, E))
.combine_data.frame <- function(x, cols) {
    if (!length(x))
        stop("length of 'x' must be > 0")
    if (!all(unlist(lapply(x, is.data.frame))))
        stop("all elements in 'x' need to be a data.frame")
    fd <- do.call(rbind, x)
    if (missing(cols))
        cols <- colnames(fd)
    else {
        if (!all(cols %in% colnames(fd)))
            stop("All columns specified with 'cols' have to be present in the",
                 " data.frames")
    }
    ## Make first unique elements
    fd <- unique(fd[, cols, drop = FALSE])
    if (nrow(fd)) {
        fd_id <- do.call(paste, fd)
        ## Find then the number of times each should be repeated.
        rep_each <- matrix(ncol = length(x), nrow = nrow(fd))
        for (i in seq_along(x)) {
            x_id <- do.call(paste, x[[i]][, cols, drop = FALSE])
            cnts <- table(match(x_id, fd_id))
            idx <- as.numeric(names(cnts))
            rep_each[as.numeric(names(cnts)), i] <- as.numeric(cnts)
        }
        splt <- split.data.frame(fd, 1:nrow(fd))
        fd <- do.call(rbind,
                      rep(splt, apply(rep_each, MARGIN = 1, max, na.rm = TRUE)))
        fd <- fd[do.call(order, fd), ]
        rownames(fd) <- NULL
    }
    fd
}

#' Convert the mzML polarity information into a character (+, -, NA)
#' representing the polarity
#'
#' @noRd
#'
#' @md
#'
#' @author Johannes Rainer
.polarity_char <- function(x) {
    x[x < 0] <- NA
    ifelse(x == 1, "+", "-")
}

#' Build an ID for a chromatogram based on the provided feature data
.chromatogramIdFromFeatureData <- function(x) {
}

readChromData <- function(files, pdata = NULL, verbose = isMSnbaseVerbose(),
                          mode = c("inMem", "onDisk")) {
    mode = match.arg(mode)
    files <- normalizePath(files)
    ## Can not use the Chromatograms object! the data is supposedly not in
    ## matrix order.
    for (i in seq_along(files)) {
        msf <- .openMSfile(files[i])
        ## This will only work for pwiz backend.
        if (!is(msf, "mzRpwiz"))
            stop("Can only extract chromatogram information from a mzML file",
                 " using the 'proteowizard' backend")
        hdr <- chromatogramHeader(msf)
        chrs <- mapply(chromatogram(msf), split.data.frame(hdr, 1:nrow(hdr)),
                       FUN = function(x, y) {
                           Chromatogram(
                               rtime = x[, 1], intensity = x[, 2],
                               precursorMz = y$precursorIsolationWindowTargetMZ,
                               productMz = y$productIsolationWindowTargetMZ,
                               fromFile = i)
                       })
        close(msf)
    }
    ## Might need an object similar to pSet:
    ## featureData that is the chromatogram info.
    ## assayData that is a list/env of Chromatogram objects (or list of lists?).
    ## There's a lot of stuff to implement!
    ## chromatogram to return a list of Chromatograms.
    ## fileIdx
    ## ...
}

## ChromExp: similar to MSnExp, just with inherited onDisk/inMem mode.

## Either extend pSet. Problem there: inherits many methods that are related to
## spectra.
## If assayData is empty but featureData is not -> onDisk mode.
## setClass("ChromExp",
##          slots = c(phenoData = "NAnnotatedDataFrame",
##                    assayData = "list",
##                    featureData = "AnnotatedDataFrame"),
##          prototype = prototype(
##              phenoData = new("NAnnotatedDataFrame"),
##              featureData = new("AnnotatedDataFrame",
##                                dimLabels = c("featureNames", "featureColumns")),
##              assayData = list()
##          )
## )

## For an object with both modes.
## - one object.
## For two objects, inMem and onDisk.
## - separate methods to access the data, otherwise I have to have if statements
## ChromExp and an OnDiskChromExp that extends ChromExp. ChromExp will contain
## THE SAME metadata, thus it will be easy to switch back and forth between them
## setClass("ChromExp",
##          contains = "pSet",
##          prototype = prototype(new("VersionedBiobase",
##                                    versions = c(classVersion("pSet"),
##                                                 ChromExp = "0.0.1"))))

## .validChromExp <- function(object) {
##     ## If it's onDisk mode: files have to be present.
##     ## If it's inMem -> check that the content of assayData fits featureData
##     ## and phenoData.
##     msg <- character()
##     ## Skip tests if assayData is empty
##     if (length(object)) {
##         ## check that nrow feature data matches assayData.
##         ## check that assayData contains only Chromatogram objects
##         ## check also validate from pSet.
##     }
    
##     if (nrow(object@featureData)) {
##     } else {
##         ## shouldn't have anything in assayData.
##     }
##     if (length(msg))
##         msg
##     else TRUE
## }

## methods to keep:
## validity???
## [
## [[
## dim
## sampleNames
## fileNames
## featureNames
## phenoData
## pData
## varMetadata
## varLabels
## featureData
## featureData<-
## fData
## fData<-
## fvarMetadata
## fvarLabels
## experimentData
## msInfo
## expinfo
## exptitle
## expemail
## ionSource
## ionSourceDetails
## analyser
## analyzer
## analyserDetails
## analyzerDetails
## instrumentModel
## instrumentManufacturer
## instrumentCustomizations
## detectorType
## description
## notes
## pubMedIds
## pubMedIds<-
## abstract
## protocolData
## processingData
## $
## $<-
## pData<-

## overwrite/replace
## precursorMz overwrite
## precursorScanNum error
## tic overwrite
## ionCount error
## precursorCharge error
## precursorIntensity error
## acquisitionNum overwrite
## scanIndex error
## rtime overwrite
## centroided error?
## smoothed error?
## peaksCount error?
## msLevel error?
## collisionEnergy overwrite
## intensity overwrite
## mz overwrite
## polarity overwrite
## fromFile overwrite
## header overwrite
## length overwrite
## spectra error
## spectrapply error

## new
## isOnDisk
