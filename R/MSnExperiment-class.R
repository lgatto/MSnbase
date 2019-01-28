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
#' @slot spectraData A [S4Vectors::DataFrame-class] storing spectra metadata.
#' @slot sampleData A [S4Vectors::DataFrame-class] storing sample metadata.
#' @slot processingQueue A `list` storing [ProcessingSteps-class] objects for
#' lazy processing.
#' @slot processing A `character` storing logging information.
#' @slot files A `character` storing absolute path to source (in general .mzML)
#' files.
#'
#' @name MSnExperiment-class
#' @docType class
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @examples
#' ## TODO
setClass(
    "MSnExperiment",
    slots=c(
        backend="Backend",
        ## was featureData in MSnExp
        spectraData="DataFrame",
        ## was phenoData in MSnExp
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
})

#' Create MSnExperiment objects
#'
#' Draft import method for the new [MSnExperiment-class]. Used as constructor.
#'
#' @param file Path to mass spectrometry data files
#' @param sampleData A [S4Vectors::DataFrame-class] object with additional
#' information to each sample (samples as rows, information as columns).
#' @param backend A [Backend-class] derivate used for internal data storage.
#' @param smoothed `logical`, are the spectra smoothed?
#' @param BPPARAM Should parallel processing be used? See
#' [BiocParallel::bpparam()].
#'
#' @rdname MSnExperiment-class
#' @export
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
    spectraData <- DataFrame(do.call(rbind, bplapply(
        file, .read_file, files=file, smoothed=smoothed, BPPARAM=BPPARAM
    )))
    backend <- backendInitialize(backend, file, spectraData, BPPARAM=BPPARAM)
    backend <- backendImportData(backend, file, spectraData, BPPARAM=BPPARAM)
    new("MSnExperiment",
        backend = backend,
        sampleData = sampleData,
        spectraData = spectraData,
        processing = paste0("Data loaded [", date(), "]"),
        files = file
    )
}
