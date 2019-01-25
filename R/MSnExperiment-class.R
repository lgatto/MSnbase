#' @include hidden_aliases.R
NULL

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
#' @name MSnExperiment-class
#' @docType class
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @examples
#' ## TODO
setClass(
    "MSnExperiment",
    slots=c(
        backend="Backend",
        ## TODO: maybe we should rename this slot to spectraData?
        featureData="DataFrame",
        ## TODO: maybe we should rename this slot to sampleData?
        ## (in SummarizedExperiment it is named colData)
        phenoData="DataFrame",
        ## Collecting ProcessingSteps for lazy processing.
        processingQueue="list",
        ## logging
        processing="character",
        files="character"
    )
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
