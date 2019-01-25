#' @include hidden_aliases.R
NULL

#' validation function for MSnExperiment
#'
#' @author Johannes Rainer
#'
#' @noRd
.validMSnExperiment <- function(object) {
    msg <- NULL
    if (length(object@files)) {
        if (length(object@files) != nrow(object@sampleData))
            msg <- c(msg, paste0("Number of files does not match the number",
                                 " of rows of 'sampleNames'"))
        if (is.null(object@backend))
            msg <- c(msg, paste0("If files are present the backend can not be",
                                 " NULL"))
    }
    if (length(msg))
        msg
    else TRUE
}

#' The MSnExperiment class
#'
#' The [MSnExperiment-class] encapsulates data and meta-data for mass
#' spectrometry experiments.
#'
#' In contrast to the old [MSnExp-class] this class supports multiple data
#' backends, e.g. in-memory ([BackendMemory-class]), on-disk as
#' mzML ([BackendMzMl-class]) or HDF5 ([BackendHdf5-class]). It supersedes
#' [MSnExp-class] and [OnDiskMSnExp-class].
#'
#' @slot backend A derivate of [Backend-class] holding/controlling the spectra
#' data.
#' @slot featureData A [S4Vectors::DataFrame-class] storing spectra metadata.
#' @slot phenoData A [S4Vectors::DataFrame-class] storing sample metadata.
#' @slot processinQueue A `list` storing [ProcessingSteps-class] objects for
#' lazy processing.
#' @slot processing A `character` storing logging information.
#' @slot files A `character` storing absolute path to source (in general .mzML)
#' files.
#'
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#'
#' @name MSnExperiment
#'
#' @examples
#' ## TODO
setClass(
    "MSnExperiment",
    slots=c(
        backend="Backend",
        spectraData="DataFrame",
        sampleData="DataFrame",
        ## Collecting ProcessingSteps for lazy processing.
        processingQueue="list",
        ## logging
        processing="character",
        files="character"
    ),
    validity = .validMSnExperiment
)


#' @rdname hidden_aliases
#' @param object Object to display.
#' @export
setMethod(
    "show",
    signature="MSnExperiment",
    definition=function(object) {
        cat("MSn experiment data (", class(object)[1L], ")\n", sep="")
        show(object@backend)
        cat("Processing:\n", paste(object@processing, collapse="\n"), "\n")
        cat("Loaded from:\n",
            paste(basename(object@files), collapse="\n"),
            "\n"
        )
})

#' @rdname MSnExperiment
readMSnExperiment <- function(file, sampleData, backend = BackendMzR(),
                              smoothed = NA, BPPARAM = bpparam()) {
    ## if (missing(backend) || !inherits(backend))
    if (missing(file) || length(file) == 0)
        stop("Parameter 'file' is required")
    if (!all(file.exists(file)))
        stop("Input file(s) can not be found")
    file <- normalizePath(file)
    if (!missing(sampleData)) {
        if (is.data.frame(sampleData))
            sampleData <- DataFrame(sampleData)
    } else {
        sampleData <- DataFrame(sampleIdx = seq_along(file))
    }
    if (!is.logical(smoothed))
        stop("smoothed should be a logical")
    .read_file <- function(z, files, smoothed) {
        file_number <- match(z, files)
        suppressPackageStartupMessages(
            require("MSnbase", quietly = TRUE, character.only = TRUE))
        msd <- .openMSfile(z)
        on.exit(close(msd))
        hdr <- header(msd)
        sp_idx <- seq_len(nrow(hdr))
        rownames(hdr) <- formatFileSpectrumNames(fileIds = file_number,
                                                 spectrumIds = seq_along(sp_idx),
                                                 nSpectra = length(sp_idx),
                                                 nFiles = length(files))
        ## rename totIonCurrent and peaksCount, as detailed in
        ## https://github.com/lgatto/MSnbase/issues/105#issuecomment-229503816
        names(hdr) <- sub("peaksCount", "originalPeaksCount", names(hdr))
        ## Add also:
        ## o fileIdx -> links to fileNames property
        ## o spIdx -> the index of the spectrum in the file.
        hdr$fileIdx <- file_number
        hdr$spIdx <- sp_idx
        hdr$smoothed <- smoothed
        if (isCdfFile(z)) {
            if (!any(colnames(hdr) == "polarity"))
                hdr$polarity <- NA
        }
        ## Order the fdData by acquisitionNum to force use of acquisitionNum
        ## as unique ID for the spectrum (issue #103). That way we can use
        ## the spIdx (is the index of the spectrum within the file) for
        ## subsetting and extracting.
        if (!all(sort(hdr$acquisitionNum) == hdr$acquisitionNum))
            warning(paste("Unexpected acquisition number order detected.",
                          "Please contact the maintainers or open an issue",
                          "on https://github.com/lgatto/MSnbase.",
                          sep = "\n")) ## see issue #160
        hdr[order(hdr$acquisitionNum), ]
    }
    new("MSnExperiment", processing = paste0("Data loaded [", date(), "]"),
        sampleData = sampleData, files = file, backend = backend,
        spectraData = DataFrame(do.call(rbind, bplapply(file, .read_file,
                                                        files = file,
                                                        smoothed = smoothed,
                                                        BPPARAM = BPPARAM))))
}
