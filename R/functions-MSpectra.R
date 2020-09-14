#' @rdname MSpectra
#'
#' @section Constructor:
#'
#' New [MSpectra] can be created with the `MSpectra(...)` function
#' where `...` can either be a single [Spectrum-class] object or a `list` of
#' `Spectrum` objects ([Spectrum1-class] and/or [Spectrum2-class]).
#'
#' @param ... For `MSpectra`: [Spectrum-class] object(s) or a `list` of
#'     [Spectrum-class] objects.
#'     For all other methods optional arguments passed along.
#'
#' @param elementMetadata For `MSpectra`: [DataFrame] with optional information
#'     that should be added as metadata information (`mcols`) to the object.
#'     The number of rows has to match the number of [Spectrum-class] objects,
#'     each row is expected to represent additional metadata information for
#'     one spectrum.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Create from Spectrum objects
#' sp1 <- new("Spectrum1", mz = c(1, 2, 4), intensity = c(4, 5, 2))
#' sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5),
#'     precursorMz = 2)
#'
#' spl <- MSpectra(sp1, sp2)
#' spl
#' spl[[1]]
#'
#' ## Add also metadata columns
#' mcols(spl)$id <- c("a", "b")
#' mcols(spl)
#'
#' ## Create a MSpectra with metadata
#' spl <- MSpectra(sp1, sp2, elementMetadata = DataFrame(id = c("a", "b")))
#'
#' mcols(spl)
#' mcols(spl)$id
MSpectra <- function(..., elementMetadata = NULL) {
    args <- list(...)
    if (length(args) == 1L && is.list(args[[1L]]))
        args <- args[[1L]]
    if (length(dim(elementMetadata)) > 1)
        rownames(elementMetadata) <- NULL
    if (is.null(names(args)))
        names(args) <- seq_len(length(args))
    new("MSpectra", listData = args, elementMetadata = elementMetadata)
}

.show_MSpectra <- function(x, margin = "", print.classinfo = FALSE) {
    cat("MSpectra with", length(x), "spectra and", length(mcols(x)),
        "metadata column(s):\n")
    if (length(x)) {
        out <- S4Vectors:::makePrettyMatrixForCompactPrinting(
                               x, .make_naked_matrix_from_MSpectra)
        if (print.classinfo) {
            .COL2CLASS <- c(msLevel = "integer", rtime = "numeric",
                            peaksCount = "integer")
            classinfo <- S4Vectors:::makeClassinfoRowForCompactPrinting(
                                         x, .COL2CLASS)
            out <- rbind(classinfo, out)
        }
        if (nrow(out))
            rownames(out) <- paste0(margin, rownames(out))
        print(out, quote = FALSE, right = TRUE, max = length(out))
    }
}

.short_spectrum_info <- function(x) {
    c(msLevel = x@msLevel, rtime = ifelse(length(x@rt), x@rt, NA_real_),
      peaksCount = peaksCount(x))
}

.make_naked_matrix_from_MSpectra <- function(x) {
    x_len <- length(x)
    mcls <- mcols(x, use.names = FALSE)
    x_mcls_len <- if (is.null(mcls)) 0L else ncol(mcls)
    res <- do.call(rbind, lapply(x, .short_spectrum_info))
    res <- apply(res, 2, format, digits = 6)
    if (!is.matrix(res))
        res <- t(as.matrix(res))
    res <- res[, c("msLevel", "rtime", "peaksCount"), drop = FALSE]
    if (x_mcls_len > 0) {
        tmp <- do.call(data.frame, c(lapply(mcls, showAsCell),
                                     list(check.names = FALSE)))
        res <- cbind(res, `|` = rep.int("|", x_len), as.matrix(tmp))
    }
    res
}

#' @title Extract data from MSnbase objects for use in Spectra
#'
#' @description
#'
#' `extractSpectraData` extracts the spectra data (m/z and intensity values
#' including metadata) from [MSnExp-class], [OnDiskMSnExp-class],
#' [Spectrum1-class], [Spectrum2-class] objects (or `list` of such objects) and
#' returns these as a `DataFrame` that can be used to create a
#' [Spectra::Spectra-class] object.This function enables thus
#' to convert data from the *old* `MSnbase` package to the newer `Spectra`
#' package.
#'
#' @param x a `list` of [Spectrum-class] objects or an object extending
#'     [MSnExp-class] or a [MSpectra-class] object.
#'
#' @return [DataFrame()] with the full spectrum data that can be passed to the
#'     [Spectra::Spectra()] function to create a `Spectra` object.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Read an mzML file with MSnbase
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_SWATH.mzML",
#'     package = "msdata")
#' data <- filterRt(readMSData(fl, mode = "onDisk"), rt = c(1, 6))
#'
#' ## Extract the data as a DataFrame
#' res <- extractSpectraData(data)
#' res
#'
#' ## This can be used as an input for the Spectra constructor of the
#' ## Spectra package:
#' ## sps <- Spectra::Spectra(res)
#' ## sps
extractSpectraData <- function(x) {
    if (inherits(x, "MSpectra")) {
        df <- DataFrame(do.call(rbind, lapply(x, .spectrum_header)))
        df <- cbind(df, mcols(x))
        df$mz <- NumericList(lapply(x, function(z) z@mz))
        df$intensity <- NumericList(lapply(x, function(z) z@intensity))
    } else if (is(x, "list") || inherits(x, "SimpleList")) {
        df <- DataFrame(do.call(rbind, lapply(x, .spectrum_header)))
        df$mz <- NumericList(lapply(x, function(z) z@mz))
        df$intensity <- NumericList(lapply(x, function(z) z@intensity))
    } else if (inherits(x, "MSnExp")) {
        df <- DataFrame(fData(x))
        df$mz <- NumericList(mz(x))
        df$intensity <- NumericList(intensity(x))
    } else stop("'x' should be either a 'list' of 'Spectrum' objects or an ",
                "object extending 'MSnExp' or 'MSpectra'.")
    colnames(df)[colnames(df) == "retentionTime"] <- "rtime"
    colnames(df)[colnames(df) == "fileIdx"] <- "fromFile"
    colnames(df)[colnames(df) == "seqNum"] <- "scanIndex"
    df
}
