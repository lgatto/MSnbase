#' @description `.writeMSData` saves the MS data from an `OnDiskMSnExp` or
#'     `MSnExp` object to mzML or mzXML files. The data is written to as many
#'     files as there are samples/files in the object. By default metadata such
#'     as data processing, original files etc are copied over from the original
#'     data files.
#'
#' @details `.writeMSData` re-calculates header columns `"peaksCount"`,
#'     `"totIonCurrent"`, `"basePeakMZ"` and `"basePeakIntensity"` on the data.
#'
#' @param object `OnDiskMSnExp` or `MSnExp` object.
#'
#' @param file `character` with the file name(s). Its length has to match the
#'     number of samples/files of `x`.
#'
#' @param outformat `character(1)` defining the format of the output files.
#'
#' @param verbose `logical(1)` if progress messages should be displayed.
#'
#' @param copy `logical(1)` if metadata (data processings, original file names
#'     etc) should be copied from the original files.
#'
#' @param software_processing optionally provide specific data processing steps.
#'     See documentation of the `software_processing` parameter of
#'     [mzR::writeMSData()].
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.writeMSData <- function(object, file, outformat = c("mzml", "mzxml"),
                         verbose = isMSnbaseVerbose(),
                         copy = TRUE,
                         software_processing = NULL) {
    ## Check input.
    outformat <- match.arg(outformat)
    if (missing(file))
        stop("'file' is required!")
    if (length(file) != length(fileNames(object)))
        stop("length of 'file' has to match the number of samples")
    if (verbose)
        message("Writing ", length(file), " ", outformat, " file.")
    if (any(file.exists(file))) {
        exst <- file[file.exists(file)]
        stop("File(s) ", paste0(exst, collapse = ", "), " does/do already exist!",
             " Please use a different file name.")
    }
    ## Split per file.
    x_split <- splitByFile(object, f = factor(fileNames(object)))
    ## Using mapply below - in principle we could then even switch to
    ## bpmapply to perform parallel saving of files.
    invisible(mapply(FUN = .writeSingleMSData, x_split, file,
                     MoreArgs = list(outformat = outformat,
                                     software_processing = software_processing,
                                     verbose = verbose, copy = copy),
                     SIMPLIFY = FALSE, USE.NAMES = FALSE))
}

#' @description `.writeMSData` for a single file.
#'
#' @md
#'
#' @noRd
.writeSingleMSData <- function(msData, file, outformat,
                               software_processing = NULL,
                               verbose = TRUE,
                               copy = is(msData, "OnDiskMSnExp")) {
    if (verbose)
        message("Saving file ", basename(file), "...", appendLF = FALSE)
    file <- path.expand(file)
    ## o get all peak data. Calling spectrapply will apply also all lazy
    ##   processing steps. We're not calling the `as` method to avoid R
    ##   having to look through all namespaces to find the function.
    pks <- spectrapply(msData, FUN = function(sp) cbind(sp@mz, sp@intensity),
                       BPPARAM = SerialParam())
        if (is(msData, "OnDiskMSnExp")) {
            hdr <- fData(msData)
            ## o Re-calculate stuff we don't have or which might have been
            ##   changed:
            ##   Add also lowMZ and highMZ; these are missing in CDF
            ##   files (issue #250).
            new_vals <- do.call(rbind, lapply(pks, function(sp) {
                max_pos <- base::which.max(sp[, 2])[1]
                pk_count <- nrow(sp)
                cbind(peaksCount = pk_count,
                      totIonCurrent = sum(sp[, 2], na.rm = TRUE),
                      basePeakMZ = sp[max_pos, 1][1],
                      basePeakIntensity = sp[max_pos, 2],
                      lowMZ = sp[1, 1],
                      highMZ = sp[pk_count, 1])
            }))
            hdr$peaksCount <- new_vals[, "peaksCount"]
            hdr$totIonCurrent <- new_vals[, "totIonCurrent"]
            hdr$basePeakMZ <- new_vals[, "basePeakMZ"]
            hdr$basePeakIntensity <- new_vals[, "basePeakIntensity"]
            hdr$lowMZ <- new_vals[, "lowMZ"]
            hdr$highMZ <- new_vals[, "highMZ"]
        } else {
            ## MSnExp: feature data does not provide all the data we need.
            hdr <- do.call(rbind, spectrapply(msData, FUN = .spectrum_header,
                                           BPPARAM = SerialParam()))
            ## In case we do have some/all fData columns: use them (issue #383).
            have_cols <- intersect(colnames(hdr), colnames(fData(msData)))
            if (length(have_cols))
                hdr[, have_cols] <- fData(msData)[, have_cols]
        }
    ## seqNum is expected to be sequentially numbered
    hdr$seqNum <- 1:nrow(hdr)
    ## Add filterString header column, if not present; export of CDF files to
    ## mzML files will otherwise fail.
    if (!any(colnames(hdr) == "filterString"))
        hdr$filterString <- NA_character_
    ## centroided has to be logical
    hdr$centroided <- as.logical(hdr$centroided)
    ## o add processing steps.
    soft_proc <- .guessSoftwareProcessing(msData, software_processing)
    if (copy) {
        copyWriteMSData(object = pks, file = file,
                        original_file = fileNames(msData), header = hdr,
                        outformat = outformat, software_processing = soft_proc)
    } else {
        writeMSData(object = pks, file = file, header = hdr,
                    outformat = outformat, software_processing = soft_proc)
    }
    if (verbose)
        message("OK")
}

#' @description Determine data processing steps based on an `OnDiskMSnExp` or
#'     `MSnExp` object and return it in a format that can be passed to the
#'     `writeMSData` or `copyWriteMSData` functions from the `mzR` package.
#'
#' @param x `OnDiskMSnExp` or `MSnExp` object.
#'
#' @param software_processing If provided it will be appended to the software
#'     processing guessed from the object. See [mzR::writeMSData()] for the
#'     format.
#'
#' @return `list` with the software processing(s).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.guessSoftwareProcessing <- function(x, software_processing = NULL) {
    res <- list()
    ## issue #321: replace MS:-1 with the MS CV Term of MSnbase, once it's
    ## included.
    msnbase_proc <- c("MSnbase", paste0(packageVersion("MSnbase"),
                                        collapse = "."), "MS:1002870")
    processings <- processingData(x)@processing
    if (length(processings)) {
        proc_cv <- unlist(lapply(processings, FUN = .pattern_to_cv))
        ## Remove all unknown processing steps (issue #321).
        proc_cv <- proc_cv[!is.na(proc_cv)]
        msnbase_proc <- c(msnbase_proc, proc_cv)
    }
    res[[1]] <- msnbase_proc
    if (length(software_processing)) {
        if (!is.list(software_processing))
            software_processing <- list(software_processing)
        res <- c(res, software_processing)
    }
    res
}


#' @description `.pattern_to_cv` performs a mapping of pattern to PSI-MS terms.
#'
#' @details The mapping bases on a manually curated list of pattern-CV-term
#'    pairs.
#'
#' @param pattern `character(1)` with the pattern for which a CV term should be
#'     returned.
#'
#' @param ifnotfound value to be returned if none of the cv terms
#'     matches the pattern.
#'
#' @param return `character` with the PSI-MS term or the value of `ifnotfound`
#'     if it was not found. Note that the length of the character can be > 1 if
#'     multiple terms match.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.pattern_to_cv <- function(pattern, ifnotfound = NA_character_) {
    if (length(pattern) > 1) {
        warning("length of pattern is > 1, using only first element")
        pattern <- pattern[1]
    }
    ## character vector mapping pattern to PSI-MS term.
    .PATTERN.TO.CV <- c(
        filter = "MS:1001486",
        normali = "MS:1001484",
        calibration = "MS:1001485",
        pick = "MS:1000035",
        centroid = "MS:1000035",
        smooth = "MS:1000592",
        baseline = "MS:1000593",
        alignment = "MS:1000745"
    )
    ## Note: we are matching each names(.PATTERN.TO.CV) to `pattern`, not the
    ## other way round. So we can provide a string describing the processing
    ## step and return the best matching CV terms.
    matches <- vapply(names(.PATTERN.TO.CV), FUN = function(z) {
        any(grepl(pattern = z, x = pattern, ignore.case = TRUE))
    }, FUN.VALUE = logical(1))
    if (any(matches))
        unname(unique(.PATTERN.TO.CV[matches]))
    else
        ifnotfound
}
