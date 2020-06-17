## Read chromatogram data from an mzML.

#' @title Read SRM/MRM chromatographic data
#'
#' @description
#'
#' The `readSRMData` function reads MRM/SRM data from provided *mzML* files and
#' returns the results as a [MChromatograms()] object.
#'
#' @details
#'
#' `readSRMData` supports reading chromatogram entries from *mzML* files. If
#' multiple files are provided the same precursor and product m/z for SRM/MRM
#' chromatograms are expected across files. The number of columns of the
#' resulting [MChromatograms()] object corresponds to the number of files. Each
#' row in the `MChromatograms` object is supposed to contain chromatograms
#' with same polarity, precursor and product m/z. If chromatograms with
#' redundant polarity, precursor and product m/z values and precursor collision
#' energies are found, they are placed into multiple consecutive rows in the
#' `MChromatograms` object.
#'
#' @note
#'
#' `readSRMData` reads only SRM/MRM chromatogram data, i.e. chromatogram data
#' from mzML files with `precursorIsolationWindowTargetMZ` and
#' `productIsolationWindowTargetMZ` attributes. Total ion chromatogram data is
#' hence not extracted.
#'
#' The number of features and hence rows of the resulting `MChromatograms`
#' object depends on the total list of unique precursor and product m/z
#' isolation windows (and precursor collision energies) found across all input
#' files. In cases in which not each file has chromatgraphic data for the same
#' polarity, precursor m/z, product m/z and collision energy,
#' an empty `Chromatogram()` object is reported for the specific precursor
#' and product m/z combination of the respective file (and a warning is
#' thrown).
#'
#' @param files `character` with the files containing the SRM/MRM data.
#'
#' @param pdata `data.frame` or `AnnotatedDataFrame` with file/sample
#'     descriptions.
#'
#' @return A `MChromatograms()` object. See details above for more information.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Read an example MRM/SRM data
#' library(msdata)
#' fl <- proteomics(full.names = TRUE, pattern = "MRM")
#'
#' ## Read the data
#' mrm <- readSRMData(fl)
#'
#' ## The data is represented as a MChromatograms object, each column
#' ## containing the data from one input file
#' mrm
#'
#' ## Access the polarity for each chromatogram (row)
#' polarity(mrm)
#'
#' ## Access the precursor m/z. The result is returned as a matrix with
#' ## columns representing the minimum and maximum m/z (will be identical in
#' ## most cases).
#' precursorMz(mrm)
#'
#' ## Access the product m/z.
#' productMz(mrm)
#'
#' ## Plot one chromatogram
#' plot(mrm[1, ])
readSRMData <- function(files, pdata = NULL) {
    files <- normalizePath(files)
    ## Read first the header of all files.
    hdr_list <- lapply(files, function(x) {
        msf <- .openMSfile(x)
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
    ## Now combine the individual headers into a feature data data.frame.
    ## What makes all this necessary is that mzML files can contain multiple
    ## chromatogram entries with exactly the same polarity, precursor and
    ## product m/z. We need to account for that.
    fdata <-  .combine_data.frame(hdr_list,
                                  cols = c("polarity",
                                           "precursorIsolationWindowTargetMZ",
                                           "productIsolationWindowTargetMZ",
                                           "precursorCollisionEnergy"))
    fdata_ids <- paste0(.polarity_char(fdata$polarity),
                       " Q1=", fdata$precursorIsolationWindowTargetMZ,
                       " Q3=", fdata$productIsolationWindowTargetMZ,
                       " collisionEnergy=", fdata$precursorCollisionEnergy)
    ## Process the phenoData.
    if (is.null(pdata))
        pdata <- data.frame(file = files, stringsAsFactors = FALSE)
    if (is.data.frame(pdata))
        pdata <- AnnotatedDataFrame(pdata)
    if (!is(pdata, "AnnotatedDataFrame"))
        stop("'pdata' should be a 'data.frame' or an 'AnnotatedDataFrame'")
    if (length(files) != nrow(pdata))
        stop("'nrow(pdata)' has to match 'length(files)'")
    pdata$file <- files

    ## Read the data.
    chrs <- unlist(mapply(
        files, hdr_list, seq_along(files),
        FUN = function(file, hdr, idx) {
            current_ids <- paste0(
                .polarity_char(hdr$polarity),
                " Q1=", hdr$precursorIsolationWindowTargetMZ,
                " Q3=", hdr$productIsolationWindowTargetMZ,
                " collisionEnergy=", hdr$precursorCollisionEnergy)
            ## Warn if the current file does contain redundant isolation windows
            if (length(current_ids) != length(unique(current_ids)))
                warning("file ", basename(file), " contains multiple ",
                        "chromatograms with identical polarity, precursor ",
                        "and product m/z values", call. = FALSE)
            ## Create the result list with the expected length
            res_chrs <- replicate(nrow(fdata), Chromatogram(fromFile = idx))
            ## Read the data
            msf <- .openMSfile(file)
            chr_data <- chromatogram(msf, hdr$chromatogramIndex)
            close(msf)
            ## Loop through all of em:
            ## 1) create Chromatogram
            ## 2) find the first emtpy slot matching the current id in res_chrs
            for (i in seq_len(nrow(hdr))) {
                idx_to_place <- which(lengths(res_chrs) == 0 &
                                      fdata_ids == current_ids[i])[1]
                if (is.na(idx_to_place))
                    stop("Got more redundant chromatograms than expected")
                res_chrs[[idx_to_place]] <- Chromatogram(
                    rtime = chr_data[[i]][, 1],
                    intensity = chr_data[[i]][, 2],
                    precursorMz = hdr$precursorIsolationWindowTargetMZ[i],
                    productMz = hdr$productIsolationWindowTargetMZ[i],
                    fromFile = idx)
            }
            res_chrs
        }), use.names = FALSE)

    MChromatograms(chrs, phenoData = pdata, featureData = fdata,
                  ncol = length(files))
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
#'
#' @noRd
#'
#' @md
#'
#' @author Johannes Rainer, Sebastian Gibb
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

    nr <- vapply(x, nrow, 1)
    dfIds <- rep(seq_along(x), nr)
    colIds <- do.call(paste, fd[, cols, drop = FALSE])

    repl <- vapply(split(dfIds, colIds), function(d)max(table(d)), 1)

    o <- order(colIds)
    fd <- fd[o,]
    colIds <- colIds[o]

    fd <- fd[rep(which(!duplicated(colIds)), repl), ]
    rownames(fd) <- NULL
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
    if (!all(x %in% c(-1, 0, 1)))
        stop("Polarity is expected to take only values 1, -1 and 0")
    x[x < 0] <- NA
    ifelse(x == 1, "+", "-")
}
