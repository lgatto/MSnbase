#' @include hidden_aliases.R
NULL

#' @title The MSnExperiment class to manage and access MS data
#'
#' @aliases MSnExperiment-class
#'
#' @name MSnExperiment
#'
#' @description
#'
#' The `MSnExperiment` class encapsules data and meta-data for mass
#' spectrometry experiments.
#'
#' In contrast to the old [MSnExp-class] this class supports multiple data
#' backends, e.g. in-memory ([BackendMemory-class]), on-disk as
#' mzML ([BackendMzMl-class]) or HDF5 ([BackendHdf5-class]). It supersedes
#' [MSnExp-class] and [OnDiskMSnExp-class] objects.
#'
#' @param all for `clean`: `logical(1)` whether all 0 intensity peaks should be
#'     removed (`TRUE`) or whether 0-intensity peaks directly adjacent to a
#'     non-zero intensity peak should be kept (`FALSE`).
#'
#' @param backend A [Backend-class] derivate used for internal data storage.
#'
#' @param BPPARAM Should parallel processing be used? See
#'     [BiocParallel::bpparam()].
#'
#' @param file `character` with the file names of the experiment.
#'
#' @param FUN for `spectrapply`: a function or the name of a function to apply
#'     to each [Spectrum-class] of the experiment.
#'
#' @param msLevel. `integer` defining the MS level of the spectra to which the
#'     function should be applied.
#'
#' @param object A `MSnExperiment` object.
#'
#' @param sampleData A [S4Vectors::DataFrame-class] object with additional
#'     information on each sample (samples as rows, information as columns).
#'
#' @param smoothed `logical`, are the spectra smoothed?
#'
#' @param t for `removePeaks`: a `numeric(1)` defining the threshold or `"min"`.
#'
#' @param verbose `logical(1)` defining the verbosity.
#'
#' @param ... for `spectrapply`: additional arguments to be passed to `FUN`.
#'
#' @section Creation of objects:
#'
#' `MSnExperiment` classes are usually created with the `readMSnExperiment`
#' function that reads general spectrum metadata information from the  mass
#' spectrometry data files.
#'
#' @section Accessing data:
#'
#' - `spectra`: get the `list` of [Spectrum-class] objects from the experiment.
#'   Note that the spectra in the `list` are not grouped by sample/file. Use
#'   the `fromFile` method to split/group the `list` by file.
#'
#' - `spectrapply`: apply an arbitrary function to each spectrum in the dataset
#'   and return its result. The function returns a `list` with the same length
#'   than there are spectra.
#'
#' @section Subsetting and filtering:
#'
#' @section Data manipulation methods:
#'
#' - `clean`: remove 0-intensity data points. See [clean()] for
#'    [Spectrum-class] objects for more details.
#'
#' - `removePeaks`: remove peaks lower than a threshold `t`. See
#'   [removePeaks()] for [Spectrum-class] objects for more details.
#'
#' @return See individual method description for the return value.
#'
#' @author Sebastian Gibb, Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## TODO
NULL

#' validation function for MSnExperiment
#'
#' @author Johannes Rainer
#'
#' @noRd
.validMSnExperiment <- function(object) {
    msg <- NULL
    if (nrow(object@sampleData)) {
        if (is.null(object@backend))
            msg <- c(msg, paste0("If sampleData is present the backend can not",
                                 " be NULL"))
        if (length(fileNames(object@backend)) != nrow(object@sampleData))
            msg <- c(msg, paste0("Number of files does not match the number",
                                 " of rows of 'sampleNames'"))
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
#'
#' @slot backend A derivate of [Backend-class] holding/controlling the spectra
#' data.
#' @slot spectraData A [S4Vectors::DataFrame-class] storing spectra metadata.
#' @slot sampleData A [S4Vectors::DataFrame-class] storing sample metadata.
#' lazy processing.
#' @slot processing A `character` storing logging information.
#'
#' @name MSnExperiment-class
#' @docType class
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setClass(
    "MSnExperiment",
    slots=c(
        backend="Backend",
        ## was featureData in MSnExp
        spectraData="DataFrame",
        ## was phenoData in MSnExp
        sampleData="DataFrame",
        ## logging
        processing="character"
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
    spectraData <- DataFrame(do.call(rbind, bplapply(
        file, .read_file, files=file, smoothed=smoothed, BPPARAM=BPPARAM
    )))
    backend <- backendInitialize(backend, file, spectraData, BPPARAM=BPPARAM)
    backend <- backendImportData(backend, spectraData, BPPARAM=BPPARAM)
    new("MSnExperiment",
        backend = backend,
        sampleData = sampleData,
        spectraData = spectraData,
        processing = paste0("Data loaded [", date(), "]")
    )
}

#' @rdname MSnExperiment
setMethod("spectrapply", "MSnExperiment", function(object, FUN = NULL,
                                                   BPPARAM = bpparam(), ...) {
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)
    isOK <- validateFeatureDataForOnDiskMSnExp(object@spectraData)
    if (length(isOK))
        stop(isOK)
    backendSpectrapply(object = object@backend,
                       spectraData = object@spectraData, FUN = FUN,
                       BPPARAM = BPPARAM, ...)
})

#' @rdname MSnExperiment
setMethod("spectra", "MSnExperiment", function(object, BPPARAM = bpparam())
    spectrapply(object = object, BPPARAM = BPPARAM))

##============================================================
##  --  DATA ACCESSORS
##
##------------------------------------------------------------

##============================================================
##  --  SUBSETTING AND FILTERING METHODS
##
##------------------------------------------------------------

##============================================================
##  --  DATA MANIPULATION METHODS
##
##------------------------------------------------------------

#' @rdname MSnExperiment
setMethod("removePeaks", "MSnExperiment", function(object, t = "min",
                                                   verbose = isMSnbaseVerbose(),
                                                   msLevel.) {
    if (!is.numeric(t) & t != "min")
        stop("Argument 't' has to be either numeric of 'min'.")
    if (missing(msLevel.))
        msLevel. <- base::sort(unique(object@spectraData$msLevel))
    if (!is.numeric(msLevel.))
        stop("'msLevel' must be numeric.")
    object@backend <- backendAddProcessing(
        object@backend, procStep = ProcessingStep("removePeaks",
                                                  list(t = t, msLevel. = msLevel.)))
    validObject(object)
    object@processing <- c(object@processing,
                           paste0("Signal <= ", t, " in MS level(s) ",
                                  paste0(msLevel., collapse = ", "),
                                  " set to 0 [", date(), "]"))
    object
})

#' @rdname MSnExperiment
setMethod("clean", "MSnExperiment", function(object, all = FALSE,
                                             verbose = isMSnbaseVerbose(),
                                             msLevel.) {
    if (!is.logical(all))
        stop("Argument 'all' must be logical")
    if (missing(msLevel.))
        msLevel. <- base::sort(unique(object@spectraData$msLevel))
    if (!is.numeric(msLevel.))
        stop("'msLevel' must be numeric.")
    object@backend <- backendAddProcessing(
        object@backend, procStep = ProcessingStep("clean",
                                                  list(all = all, msLevel. = msLevel.)))
    validObject(object)
    object@processing <- c(object@processing,
                           paste0("Spectra of MS level(s) ",
                                  paste0(msLevel., collapse = ", "),
                                  " cleaned [", date(), "]"))
    object
})
