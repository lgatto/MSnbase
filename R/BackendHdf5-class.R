#' @include hidden_aliases.R
NULL

#' Backend class that stores data temporarily to hdf5 files.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @md
#'
#' @noRd
setClass(
    "BackendHdf5",
    contains="Backend",
    slots=c(checksums="character", h5files="character")
)

#' @rdname hidden_aliases
setMethod("show", "BackendHdf5", function(object) {
    callNextMethod()
    cat("Source files:\n",
        paste(" ", unique(dirname(object@h5files)), collapse="\n"), "\n", sep=""
    )
})

#' @rdname Backend
BackendHdf5 <- function() new("BackendHdf5")

.valid.BackendHdf5.checksums <- function(checksums, h5files) {
    if (length(checksums) != length(h5files))
        return("different number of md5 sums and hdf5 files")
    chk <- .vdigest(h5files, file=TRUE)
    m <- checksums == chk
    if (any(!m)) {
        return(paste0(
            "Checksum(s) for file(s) don't match:\n",
            paste0(
                " ", h5files[!m],
                " (current: ", chk[!m], ", stored: ", checksums[!m], ")",
                collapse="\n"
            )
        ))
    }
    NULL
}

.valid.BackendHdf5.h5files <- function(src, h5files) {
    if (length(src) != length(h5files))
        return("Different number of hdf5 and original source files.")
    if (anyDuplicated(h5files))
        return("Hdf5 files are not unique.")
    if (!all(file.exists(h5files))) {
        return(paste0(
            "File(s) ",
            paste0(h5files[!file.exists(h5files)], collapse=", "),
            " not found."
        ))
    }
    NULL
}

setValidity("BackendHdf5", function(object) {
    msg <- c(
        .valid.BackendHdf5.h5files(object@files, object@h5files),
        .valid.BackendHdf5.checksums(object@checksums, object@h5files)
    )
    if (is.null(msg)) { TRUE } else { msg }
})

#' Initialize the Hdf5 backend, i.e. create the hdf5 files and add their names
#' to the @hdf5file slot.
#'
#' @param object `BackendHdf5`
#'
#' @param files `character` with the original file names.
#'
#' @param spectraData `DataFrame` with the spectrum metadata.
#'
#' @param path `character(1)` defining the path where the hdf5 files should be
#'     stored.
#'
#' @return `BackendHdf5` object.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @md
#'
#' @noRd
setMethod(
    "backendInitialize",
    "BackendHdf5",
    function(object, files, spectraData, path=".", ...) {
    path <- normalizePath(path, mustWork=FALSE)
    dir.create(path, showWarnings=FALSE, recursive=TRUE)

    files <- normalizePath(files)
    object@h5files <- file.path(path, paste0(.vdigest(files), ".h5"))

    if (any(file.exists(object@h5files))) {
        stop(
            "File(s) ",
            paste0(object@h5files[file.exists(object@h5files)], collapse=", "),
            " already exist(s). Please choose a different 'path'."
        )
    }

    for (i in seq(along=object@h5files)) {
        h5 <- H5Fcreate(object@h5files[i])
        h5createGroup(h5, "spectra")
        H5Fclose(h5)
    }

    object@checksums <- .vdigest(object@h5files, file=TRUE)

    callNextMethod()
})

#' Import the data from the raw MS files and store them to the hdf5 files.
#'
#' @inheritParams backendInitialize
#'
#' @return `BackendHdf5`
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @md
#'
#' @noRd
setMethod(
    "backendImportData",
    "BackendHdf5",
    function(object, spectraData, ..., BPPARAM=bpparam()) {

    bpmapply(object@files, object@h5files, FUN=.mzr2hdf5, BPPARAM=BPPARAM)
    object@checksums <- .vdigest(object@h5files, file=TRUE)
    validObject(object)
    object
})

setMethod(
    "backendReadSpectra",
    "BackendHdf5",
    function(object, spectraData, ...) {

    uFileIdx <- unique(spectraData$fileIdx)
    h5files <- object@h5files[uFileIdx]

    changed <- .valid.BackendHdf5.checksums(object@checksums[uFileIdx], h5files)
    if (!is.null(changed))
        stop(changed)

    spectra <- vector(mode="list", length=nrow(spectraData))
    spectra <- split(spectra, spectraData$fileIdx)
    spd <- split(spectraData, spectraData$fileIdx)
    spIdx <- split(spectraData$spIdx, spectraData$fileIdx)

    for (i in seq(along=h5files)) {
        h5fh <- H5Fopen(h5files[i])
        h5gfh <- H5Gopen(h5fh, "spectra")
        spectra[[i]] <- .spectra_from_data(
            lapply(as.character(spIdx[[i]]), h5read, file=h5gfh), spd[[i]]
        )
        H5Gclose(h5gfh)
        H5Fclose(h5fh)
    }

    spectra <- unsplit(spectra, spectraData$fileIdx)
    names(spectra) <- rownames(spectraData)
    spectra
})

setMethod(
    "backendWriteSpectra",
    "BackendHdf5",
    function(object, spectra, spectraData) {

    lvl <- .hdf5_compression_level()

    uFileIdx <- unique(spectraData$fileIdx)
    h5files <- object@h5files[uFileIdx]
    spectra <- split(spectra, spectraData$fileIdx)

    for (i in seq(along=h5files)) {
        h5fh <- H5Fopen(h5files[i])
        h5gfh <- H5Gopen(h5fh, "spectra")
        for (j in seq(along=spectra[[i]])) {
            h5write(
                cbind(mz(spectra[[i]][[j]]), intensity(spectra[[i]][[j]])),
                h5gfh, as.character(i), level=lvl
            )
        }
        H5Gclose(h5gfh)
        H5Fclose(h5fh)
    }

    object@checksums[uFileIdx] <- .vdigest(h5files, file=TRUE)
    validObject(object)
    object
})

#' Write the content of a single mzML/etc file to an h5file. We're using the
#' spectrum index in the file as data set ID.
#'
#' @return `character(1)` with the md5 sum of the spectra data.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @noRd
.mzr2hdf5 <- function(file, h5) {
    msfh <- openMSfile(file)
    ## read complete header first to avoid random seq fault errors from
    ## proteowizard
    hdr <- header(msfh)
    on.exit(close(msfh))

    lvl <- .hdf5_compression_level()

    h5fh <- H5Fopen(h5)
    h5gfh <- H5Gopen(h5fh, "spectra")
    on.exit(H5Gclose(h5gfh), add=TRUE)
    on.exit(H5Fclose(h5fh), add=TRUE)

    for (i in seq(along=msfh)) {
        h5write(peaks(msfh, i), h5gfh, as.character(i), level=lvl)
    }
}
