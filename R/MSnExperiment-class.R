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
#' It supersedes [MSnExp-class] and [OnDiskMSnExp-class] objects and supports
#' multiple data backends, e.g. in-memory ([BackendMemory-class]), on-disk as
#' mzML ([BackendMzR-class]) or HDF5 ([BackendHdf5-class]).
#'
#' @details
#'
#' The `MSnExperiment` class uses by default a lazy data manipulation strategy,
#' i.e. data manipulations such as performed with `removePeaks` are not applied
#' immediately to the data, but applied on-the-fly to the spectrum data once it
#' is retrieved.
#'
#' @param all for `clean`: `logical(1)` whether all 0 intensity peaks should be
#'     removed (`TRUE`) or whether 0-intensity peaks directly adjacent to a
#'     non-zero intensity peak should be kept (`FALSE`).
#'
#' @param backend a [Backend-class] derivate used for internal data storage.
#'
#' @param BPPARAM should parallel processing be used? See
#'     [BiocParallel::bpparam()].
#'
#' @param drop for `[`: if `drop = TRUE` and the object is subsetted to a single
#'     element, a `Spectrum` class is returned; `drop = FALSE` returns always
#'     a `MSnExperiment` object.
#'
#' @param f for `spectrapply`: `factor`, `character`, `numeric` or `logical`
#'     (same length than there are spectra in `object`, i.e. with length
#'     equal to `nrow(spectraData(object))` to define how the data should be
#'     split into chunks for parallelization.
#'
#' @param file for `readMSnExperiment: `character` with the file names of the
#'     experiment. For `filterFile`: index or name of the file to which the
#'     data should be subsetted.
#'
#' @param FUN for `spectrapply`: a function or the name of a function to apply
#'     to each [Spectrum-class] of the experiment.
#'
#' @param i for `[`: `integer`, `logical` or `character` specifying the
#'     **spectra** to which `object` should be subsetted.
#'
#' @param j for `[`: not supported.
#'
#' @param metadata for `MSnExperiment` and `readMSnExperiment`: `list` with
#'     optional metadata information.
#'
#' @param msLevel. `integer` defining the MS level of the spectra to which the
#'     function should be applied.
#'
#' @param object a `MSnExperiment` object.
#'
#' @param sampleData a [S4Vectors::DataFrame-class] object with additional
#'     information on each sample (samples as rows, information as columns).
#'
#' @param spectraData for `MSnExperiment`: a [S4Vectors::DataFrame-class] object
#'     with optional additional metadata columns for each spectrum.
#'
#' @param smoothed `logical`, are the spectra smoothed?
#'
#' @param t for `removePeaks`: a `numeric(1)` defining the threshold or `"min"`.
#'
#' @param value for `featureData`, `sampleData` and `spectraData`: a `DataFrame`
#'
#' @param verbose `logical(1)` defining the verbosity.
#'
#' @param x a `MSnExperiment` object.
#'
#' @param ... for `readMSnExperiment`: additional parameters to be passed to
#'     the init method of the backend, such as `path` for [BackendHdf5()] to
#'     define the directory where the hdf5 files should be saved.
#'     For `spectrapply`: additional arguments to be passed to `FUN`.
#'
#' @section Creation of objects, conversion and changing the backend:
#'
#' `MSnExperiment` classes are usually created with the `readMSnExperiment`
#' function that reads general spectrum metadata information from the  mass
#' spectrometry data files.
#'
#' Alternatively it is possible to create a new object from a list of `Spectrum`
#' objects using the `MSnExperiment` function. Additional spectrum metadata
#' columns can be provided with the `spectraData` argument, sample annotations
#' with the `sampleData` argument and arbitrary metadata with the `metadata`
#' argument. Note that objects created with the `MSnExperiment` constructor
#' function can not use the `BackendMzR` as backend.
#'
#' `MSnExperiment` objects can be converted to a `list` or
#' [S4Vectors::List-class] of `Spectrum` objects with the `as(object, "list")`
#' and `as(object, "List")` function, respectively.
#'
#' The [Backend-class] can be changed with the `setBackend` function by
#' specifying the new [Backend-class] with the `backend` parameter. See examples
#' for more details.
#'
#' @section Accessing data:
#'
#' - `acquisitionNum`: get the acquisition number of each spectrum as a
#'   named `integer` vector with the same length than `object`.
#'
#' - `featureData`: get or set general spectrum metadata. Returns a `DataFrame`
#'   or a `MSnExperiment` with updated spectra metadata. Each row of the
#'   `DataFrame` contains information for one spectrum. This function is
#'   equivalent to [featureData()] of `MSnExp`/`OnDiskMSnExp` objects.
#'
#' - `featureNames`: extract the feature (spectrum) names.
#'
#' - `fileNames`: get the original file names from which the data was imported.
#'
#' - `length`: get the number of spectra in the object.
#'
#' - `metadata`: get the metadata `list`.
#'
#' - `sampleData`: get or set sample metadata. Returns a `DataFrame`, each row
#'   containing information for one sample or file or a `MSnExperiment` with
#'   the update sample metadata. This function is equivalent to [phenoData()]
#'   of `MSnExp`/`OnDiskMSnExp` objects.
#'
#' - `spectraData`: get or set general spectrum metadata. See `featureData`
#'   above.
#'
#' - `spectrapply`: apply an arbitrary function to each spectrum in the dataset
#'   and return its result. The function returns a `list` with the same length
#'   than there are spectra. Argument `f` allows to define how to split the
#'   data/spectra into chunks for paralellization. By default data access and
#'   application of the provided function are parallelized by file.
#'
#' @section Subsetting and filtering:
#'
#' - `[i]`: subset the object by spectra (`i`). Returns an `MSnExperiment`,
#'   unless `drop = TRUE` and the object is subsetted to a single spectrum,
#'   in which case a `Spectrum` is returned.
#'
#' - `[[i]]`: extract the [Spectrum-class] with index `i` from the data.
#'
#' - `filterFile`: subset the object by file. Returns an `MSnExperiment`.
#'
#' @section Data manipulation methods:
#'
#' - `clean`: remove 0-intensity data points. See [clean()] for
#'    [Spectrum-class] objects for more details.
#'
#' - `removePeaks`: remove peaks lower than a threshold `t`. See
#'   [removePeaks()] for [Spectrum-class] objects for more details.
#'
#'
#' @return See individual method description for the return value.
#'
#' @author Sebastian Gibb, Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Create an MSnExperiment from a list of Spectrum objects
#' spl <- list(new("Spectrum1", rt = 1.2, mz = 1:4, intensity = abs(rnorm(4))),
#'     new("Spectrum1", rt = 1.3, mz = 1:4, intensity = abs(rnorm(4))))
#' mse <- MSnExperiment(spl)
#' mse
#'
#' ## Access the second spectrum
#' mse[[2]]
#'
#' ## Get the spectrum metadata
#' spectraData(mse)
#'
#' ## Add an additional column to the spectrum metadata
#' spectraData(mse)$peak_id <- c("a", "b")
#'
#' ## featureData and spectraData both access the spectrum metadata
#' featureData(mse)
#'
#'
#' ## Create an MSnExperiment from two input files using the on-disk
#' ## BackendMzR backend
#' sf <- dir(system.file("sciex", package = "msdata"), full.names = TRUE)
#' mse <- readMSnExperiment(sf, backend = BackendMzR())
#' mse
#'
#' ## Get associated file names
#' fileNames(mse)
#'
#' ## Get spectra metadata
#' spectraData(mse)
#'
#' ## Subset the object to contain only spectra 3, 12, 45
#' mse_sub <- mse[c(3, 12, 45)]
#' mse_sub
#'
#' ## Coerce to a list of spectra
#' as(mse_sub, "list")
#'
#' ## Subset the object to contain only spectra from the second file
#' mse_sub <- filterFile(mse, 2)
#' fileNames(mse_sub)
#'
#' ## Apply an arbitrary function to each spectrum and return its results.
#' ## We calculate the mean intensity for each spectrum. By
#' ## default the function parallelizes the operation FUN by file.
#' res <- spectrapply(mse, FUN = function(z) mean(intensity(z)))
#' head(res)
#'
#' ## Parameter `f` can be used to specify how the function splits the data
#' ## into chunks for parallel processing. Below we disable parallel processing
#' ## by defining a single chunk.
#' res <- spectrapply(mse, f = rep(1, nrow(spectraData(mse))),
#'     FUN = function(z) mean(intensity(z)))
#' head(res)
#'
#' ## The `setBackend` function can be used to change the backend for the
#' ## `MSnExperiment`. Below we change the backend from the default raw MS
#' ## data files-based backend (`BackendMzR`) to the HDF5-file based
#' ## `BackendHdf5`. With the additional `path` parameter we specify the
#' ## directory in which the HDF5 files should be saved.
#' mse <- setBackend(mse, backend = BackendHdf5(),
#'     path = paste0(tempdir(), "/hdf5"))
#' mse
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
    msg <- c(msg, .valid.processingQueue(object@processingQueue),
             .valid.MSnExperiment.featureData(object@spectraData,
                                              nrow(object@sampleData)))
    if (length(msg))
        msg
    else TRUE
}

.valid.processingQueue <- function(x) {
    if (length(x))
        if (!all(vapply(x, inherits, logical(1), "ProcessingStep")))
            return("'processingQueue' should only contain ProcessingStep objects.")
    NULL
}

.valid.MSnExperiment.featureData <- function(x, nsamples) {
    msg <- NULL
    if (nrow(x)) {
        if (!any(colnames(x) == "fileIdx"))
            msg <- c(msg, "Column 'fileIdx' is required")
        else if (anyNA(x$fileIdx) | !is.integer(x$fileIdx))
            msg <- c(msg, "'fileIdx' should contain only non-missing integers")
        if (anyDuplicated(rownames(x)))
            msg <- c(msg, "rownames have to be unique")
        if (!missing(nsamples))
            if (length(unique(x$fileIdx)) != nsamples)
                msg <- c(msg, paste0("number of distinct indices in 'fileIdx' ",
                                     "does not match number of files/samples"))
        if (!any(colnames(x) == "msLevel"))
            msg <- c(msg, "Column 'msLevel' is required")
        else if (anyNA(x$msLevel) | !is.integer(x$msLevel))
            msg <- c(msg, "'msLevel' should contain only non-missing integers")
    }
    if (length(msg))
        msg
    else NULL
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
#' @slot processingQueue `list` of `ProcessingStep` objects.
#' @slot processing A `character` storing logging information.
#' @slot metadata A `list` storing experiment metadata.
#' @slot version A `characher(1)` containing the class version.
#'
#' @name MSnExperiment-class
#' @docType class
#' @author Sebastian Gibb \email{mail@@sebastiangibb.de}
#' @noRd
setClass(
    "MSnExperiment",
    slots = c(
        backend = "Backend",
        ## was featureData in MSnExp
        spectraData = "DataFrame",
        ## was phenoData in MSnExp
        sampleData = "DataFrame",
        processingQueue = "list",
        ## logging
        processing = "character",
        ## metadata
        metadata = "list",
        version = "character"
    ),
    prototype = prototype(version = "0.1"),
    validity = .validMSnExperiment
)

#' @rdname MSnExperiment
MSnExperiment <- function(x, spectraData, sampleData, metadata, ...) {
    if (missing(x) || any(vapply(x, function(z) !inherits(z, "Spectrum"),
                                 logical(1))))
        stop("'x' has to be a list of 'Spectrum' objects")
    fdata <- DataFrame(do.call(rbind, lapply(x, .spectrum_header)))
    if (!missing(spectraData))
        spectraData <- .merge_featureData(spectraData, fdata)
    else spectraData <- fdata
    if (all(is.na(spectraData$fileIdx)))
        spectraData$fileIdx <- 1L
    file <- unique(spectraData$fileIdx)
    if (!is.null(names(x)))
        rownames(spectraData) <- names(x)
    else {
        names(x) <- formatFileSpectrumNames(
            fileIds = spectraData$fileIdx,
            spectrumIds = unlist(lapply(table(
                factor(spectraData$fileIdx, unique(spectraData$fileIdx))
            ), seq_len), use.names = FALSE),
            nFiles = length(file))
        rownames(spectraData) <- names(x)
    }
    if (missing(sampleData))
        sampleData <- DataFrame(sampleIdx = unique(fdata$fileIdx))
    backend <- BackendMemory()
    backend <- backendInitialize(backend, as.character(file), spectraData, ...)
    backend <- backendWriteSpectra(backend, x, spectraData)
    new("MSnExperiment",
        backend = backend,
        sampleData = sampleData,
        spectraData = spectraData,
        processingQueue = list(),
        processing = paste0("Data loaded [", date(), "]")
        )
}

#' @description
#'
#' Merge two featureData `DataFrame`s. The resulting `DataFrame` contains all
#' columns from both `x` and `y`. For columns present in both `DataFrame`s those
#' in `x` will be used. Also, the resulting `DataFrame` uses the row names of
#' `x` unless `x` has no row names.
#'
#' @param x `DataFrame`
#'
#' @param y `DataFrame`
#'
#' @return `DataFrame` with the merged columns.
#'
#' @author Johannes Rainer
#'
#' @noRd
.merge_featureData <- function(x, y) {
    if (nrow(x) != nrow(y))
        stop("'x' and 'y' have to have the same number of rows")
    if (is.null(rownames(x)) & !is.null(rownames(y)))
        rownames(x) <- rownames(y)
    cols_y <- !(colnames(y) %in% colnames(x))
    if (any(cols_y))
        x <- cbind(x, y[, cols_y, drop = FALSE])
    x
}

#' @rdname hidden_aliases
#' @param object Object to display.
#' @export
setMethod(
    "show",
    signature="MSnExperiment",
    definition=function(object) {
        cat("MSn data (", class(object)[1L], ") with ",
            nrow(object@spectraData), " spectra:\n", sep="")
        txt <- capture.output(
            object@spectraData[, c("msLevel", "retentionTime", "totIonCurrent")])
        cat(txt[-1], sep = "\n")
        show(object@backend)
        if (length(object@processingQueue))
            cat("Lazy evaluation queue:", length(object@processingQueue),
                "processing step(s)\n")
        cat("Processing:\n", paste(object@processing, collapse="\n"), "\n")
    })

#' @rdname MSnExperiment
readMSnExperiment <- function(file, sampleData, backend = BackendMzR(),
                              smoothed = NA, metadata = list(), ...,
                              BPPARAM = bpparam()) {
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
    spectraData <- DataFrame(
        do.call(rbind, bplapply(file, .read_file, files=file,
                                smoothed=smoothed, BPPARAM=BPPARAM)))
    backend <- backendInitialize(backend, file, spectraData, ...,
                                 BPPARAM=BPPARAM)
    backend <- backendImportData(backend, spectraData, BPPARAM=BPPARAM)
    new("MSnExperiment",
        backend = backend,
        sampleData = sampleData,
        spectraData = spectraData,
        processingQueue = list(),
        metadata = metadata,
        processing = paste0("Data loaded [", date(), "]")
    )
}

#' @rdname MSnExperiment
setGeneric("setBackend", function(object, backend, ..., BPPARAM = bpparam())
    standardGeneric("setBackend"))
#' @rdname hidden_aliases
setMethod("setBackend", c("MSnExperiment", "Backend"),
          function(object, backend, ..., BPPARAM = bpparam()) {
              backend <- backendInitialize(backend, fileNames(object),
                                           object@spectraData, ...)
              backend <- backendWriteSpectra(
                  backend, backendReadSpectra(object@backend,
                                              object@spectraData),
                  object@spectraData)
              object@backend <- backend
              validObject(object)
              object
          })
#' @rdname hidden_aliases
setMethod("setBackend", c("MSnExperiment", "BackendMzR"),
          function(object, backend, ..., BPPARAM = bpparam()) {
              ## TODO: add check for change in backend once implemented by
              ## @sgibb
              object@backend <- backendInitialize(backend, fileNames(object),
                                                  object@spectraData, ...)
              validObject(object)
              object
          })
#' @rdname hidden_aliases
setMethod("setBackend", c("MSnExperiment", "BackendHdf5"),
          function(object, backend, ..., BPPARAM = bpparam()) {
              backend <- backendInitialize(backend, fileNames(object),
                                           object@spectraData, ...)
              spd <- split(object@spectraData, object@spectraData$fileIdx)
              cnts <- bplapply(spd, function(z, hdf5_backend, backend) {
                  res <- backendWriteSpectra(hdf5_backend,
                                             backendReadSpectra(backend, z), z)
                  res@checksums[z$fileIdx[1]]
              }, hdf5_backend = backend, backend = object@backend,
              BPPARAM = BPPARAM)
              backend@checksums <- unlist(cnts)
              object@backend <- backend
              validObject(object)
              object
})


#' @rdname MSnExperiment
setMethod("spectrapply", "MSnExperiment", function(object,
                                                   f = spectraData(object)$fileIdx,
                                                   FUN = NULL,
                                                   BPPARAM = bpparam(), ...) {
    BPPARAM <- getBpParam(object, BPPARAM = BPPARAM)
    isOK <- validateFeatureDataForOnDiskMSnExp(object@spectraData)
    if (length(isOK))
        stop(isOK)
    if (length(f) != nrow(spectraData(object)))
        stop("Length of 'f' has to match 'nrow(spectraData(object))'")
    pqueue <- object@processingQueue
    if (!is.null(FUN))
        pqueue <- c(pqueue, ProcessingStep(FUN, ARGS = list(...)))
    res <- bplapply(split(object@spectraData, f), function(z, queue, bknd) {
        .apply_processing_queue(backendReadSpectra(bknd, z), queue)
    }, queue = pqueue, bknd = object@backend, BPPARAM = BPPARAM)
    names(res) <- NULL
    res <- unlist(res, recursive = FALSE)
    res[rownames(object@spectraData)]
})

#' Helper function to add an arbitrary function with its arguments as a
#' processing step to the object's `processingQueue`.
#'
#' @param object any object with an `processingQueue` slot.
#'
#' @param FUN function or name of a function.
#'
#' @param ... Additional arguments to `FUN`.
#'
#' @author Johannes Rainer
#'
#' @noRd
addProcessingStep <- function(object, FUN, ...) {
    if (missing(FUN))
        return(object)
    object@processingQueue <- c(object@processingQueue,
                                list(ProcessingStep(FUN, ARGS = list(...))))
    validObject(object)
    object
}

setAs("MSnExperiment", "list", function(from) {
    spectrapply(from)
})
setAs("MSnExperiment", "List", function(from) {
    List(spectrapply(from))
})

##============================================================
##  --  DATA ACCESSORS
##
##------------------------------------------------------------

#' @rdname MSnExperiment
setMethod("featureData", "MSnExperiment", function(object) {
    object@spectraData
})

#' @rdname MSnExperiment
setReplaceMethod("featureData", "MSnExperiment", function(object, value) {
    if (!is(value, "DataFrame"))
        stop("'value' should be a 'DataFrame'")
    if (nrow(value) != nrow(object@spectraData))
        stop("Expecting ", nrow(object@spectraData), " rows, but 'value' has",
             " only ", nrow(value))
    object@spectraData <- value
    validObject(object) # check before updating backend
    object@backend <- backendUpdateMetadata(object@backend, value)
    object
})

#' @rdname MSnExperiment
setMethod("spectraData", "MSnExperiment", function(object) {
    featureData(object)
})

#' @rdname MSnExperiment
setReplaceMethod("spectraData", "MSnExperiment", function(object, value) {
    featureData(object) <- value
    object
})

#' @rdname MSnExperiment
setMethod("sampleData", "MSnExperiment", function(object) {
    object@sampleData
})

#' @rdname MSnExperiment
setReplaceMethod("sampleData", "MSnExperiment", function(object, value) {
    if (!is(value, "DataFrame"))
        stop("'value' should be a 'DataFrame'")
    object@sampleData <- value
    validObject(object)
    object
})

#' @rdname MSnExperiment
setMethod("metadata", "MSnExperiment",
          function(x, ...) {
              if (is.null(x@metadata) || is.character(x@metadata))
                  list(metadata = x@metadata)
              else x@metadata
          })

#' @rdname MSnExperiment
setMethod("fileNames", "MSnExperiment", function(object) {
    fileNames(object@backend)
})

#' @rdname MSnExperiment
setMethod("acquisitionNum", "MSnExperiment", function(object) {
    res <- if (is.null(object@spectraData$acquisitionNum))
               rep_len(NA_integer_, length(object))
           else object@spectraData$acquisitionNum
    names(res) <- featureNames(object)
    res
})

## centroided
## collisionEnergy
#' @rdname MSnExperiment
setMethod("featureNames", "MSnExperiment", function(object) {
    rownames(object@spectraData)
})
## fromFile
## intensity
## ionCount
## isCentroided
## isEmpty
#' @rdname MSnExperiment
setMethod("length", "MSnExperiment", function(x) {
    nrow(x@spectraData)
})
## msLevel
## mz
## polarity
## rtime
## peaksCount
## precursorCharge
## precursorIntensity
## precursorMz
## precScanNum
## scanIndex
## smoothed
## tic

##============================================================
##  --  SUBSETTING AND FILTERING METHODS
##
##------------------------------------------------------------
#' @rdname MSnExperiment
setMethod("[", "MSnExperiment", function(x, i, j, ..., drop = TRUE) {
    if (!missing(j))
        stop("Subsetting by columns/samples is not supported")
    if (missing(i))
        return(x)
    i <- .to_index(rownames(x@spectraData), i)
    x@spectraData <- x@spectraData[i, , drop = FALSE]
    x@backend <- backendSubset(x@backend, x@spectraData)
    file <- unique(x@spectraData$fileIdx) # here we allow unsorted file idx.
    x@spectraData$fileIdx <- match(x@spectraData$fileIdx, file)
    x@sampleData <- x@sampleData[file, , drop = FALSE]
    if (nrow(x@spectraData) == 1 & drop)
        x <- spectrapply(x)[[1]]
    validObject(x)
    x
})

#' @rdname MSnExperiment
setMethod("[[", "MSnExperiment",
          function(x, i, j = "missing", drop = "missing") {
              x[i, , drop = TRUE]
          })

#' @rdname MSnExperiment
setMethod("filterFile", "MSnExperiment", function(object, file) {
    if (missing(file))
        return(object)
    file <- .to_index(fileNames(object), file, variable = "file")
    object@spectraData <- object@spectraData[object@spectraData$fileIdx %in%
                                             file, , drop = FALSE]
    object@spectraData <- object@spectraData[order(match(
                                     object@spectraData$fileIdx, file)), ,
                                     drop = FALSE]
    object@backend <- backendSubset(object@backend, object@spectraData)
    object@spectraData$fileIdx <- match(object@spectraData$fileIdx, file)
    object@sampleData <- object@sampleData[file, , drop = FALSE]
    object@processing <- c(object@processing,
                           paste0("Filter: select file(s): ",
                                  paste0(file, collapse = ", "),
                                  " [", date(), "]"))
    validObject(object)
    object
})

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
    object <- addProcessingStep(object, "removePeaks", t = t,
                                msLevel. = msLevel.)
    object@processing <- c(object@processing,
                           paste0("Signal <= ", t, " in MS level(s) ",
                                  paste0(msLevel., collapse = ", "),
                                  " set to 0 [", date(), "]"))
    validObject(object)
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
    object <- addProcessingStep(object, "clean", all = all, msLevel. = msLevel.)
    object@processing <- c(object@processing,
                           paste0("Spectra of MS level(s) ",
                                  paste0(msLevel., collapse = ", "),
                                  " cleaned [", date(), "]"))
    validObject(object)
    object
})
