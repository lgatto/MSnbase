#' @description `.writeMSData` saves the MS data from an `OnDiskMSnExp` or
#'     `MSnExp` object to mzML or mzXML files. The data is written to as many
#'     files as there are samples/files in the object. By default metadata such
#'     as data processing, original files etc are copied over from the original
#'     data files.
#'
#' @details `.writeMSData` re-calculates header columns `"peaksCount"`,
#'     `"totIonCurrent"`, `"basePeakMZ"` and `"basePeakIntensity"` on the data.
#'
#' @param x `OnDiskMSnExp` or `MSnExp` object.
#'
#' @param files `character` with the file name(s). Its length has to match the
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
.writeMSData <- function(x, files, outformat = c("mzml", "mzxml"),
                         verbose = isMSnbaseVerbose(),
                         copy = TRUE,
                         software_processing = NULL) {
    ## Check input.
    outformat <- match.arg(outformat)
    if (missing(files))
        stop("'files' is required!")
    if (length(files) != length(fileNames(x)))
        stop("length of 'files' has to match the number of samples")
    if (verbose)
        message("Writing ", length(files), " ", outformat, " files.")
    ## Split per file.
    x_split <- splitByFile(x, f = factor(fileNames(x)))
    ## Using mapply below - in principle we could then even switch to
    ## bpmapply to perform parallel saving of files.
    res <- mapply(x_split, files, FUN = function(z, f, outformat,
                                                 software_processing,
                                                 verbose){
        if (verbose)
            message("Saving file ", basename(f), "...", appendLF = FALSE)
        if (file.exists(f))
            stop("File ", f, " does already exist! Please use a different ",
                 "file name.")
        ## o get all peak data. Calling spectrapply will apply also all lazy
        ##   processing steps. We're not calling the `as` method to avoid R
        ##   having to look through all namespaces to find the function.
        pks <- spectrapply(z, FUN = function(sp) cbind(sp@mz, sp@intensity),
                           BPPARAM = SerialParam())
        if (is(x, "OnDiskMSnExp")) {
            hdr <- fData(z)
        } else {
            ## MSnExp: feature data does not provide all the data we need.
            hdr <- data.frame(do.call(rbind,
                                      spectrapply(z, FUN = .spectrum_header,
                                                  BPPARAM = SerialParam())),
                              check.names = FALSE)
        }
        ## o Re-calculate stuff we don't have or which might have been
        ##   changed:
        new_vals <- do.call(rbind, lapply(pks, function(sp) {
            cbind(peaksCount = nrow(sp),
                  totIonCurrent = sum(sp[, 2], na.rm = TRUE),
                  basePeakMZ = sp[base::which.max(sp[, 2]), 1][1],
                  basePeakIntensity = max(sp[, 2], na.rm = TRUE))
        }))
        hdr$peaksCount <- new_vals[, "peaksCount"]
        hdr$totIonCurrent <- new_vals[, "totIonCurrent"]
        hdr$basePeakMZ <- new_vals[, "basePeakMZ"]
        hdr$basePeakIntensity <- new_vals[, "basePeakIntensity"]
        ## seqNum is expected to be sequentially numbered
        hdr$seqNum <- 1:nrow(hdr)
        ## o add processing steps.
        soft_proc <- .guessSoftwareProcessing(z, software_processing)
        if (copy) {
            copyWriteMSData(f, original_file = fileNames(z), header = hdr,
                            data = pks, outformat = outformat,
                            software_processing = soft_proc)
        } else {
            writeMSData(f, header = hdr, data = pks, outformat = outformat,
                        software_processing = soft_proc)
        }
        if (verbose)
            message("OK")
        NULL
        }, MoreArgs = list(outformat = outformat,
                           software_processing = software_processing,
                           verbose = verbose))
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
    msnbase_proc <- c("MSnbase", paste0(packageVersion("MSnbase"),
                                        collapse = "."))
    processings <- processingData(x)@processing
    if (length(processings)) {
        proc_cv <- unlist(lapply(processings, FUN = .pattern_to_cv))
        ## Remove MS:-1 if we have also some mapped processings.
        if (!all(proc_cv == "MS:-1"))
            proc_cv <- proc_cv[proc_cv != "MS:-1"]
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
#' @param ifnotfound `character(1)` to be returned if none of the cv terms
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
.pattern_to_cv <- function(pattern, ifnotfound = "MS:-1") {
    res <- ifnotfound
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
        smooth = "MS:1000542",
        baseline = "MS:1000543",
        alignment = "MS:1000745"
    )
    ## Note: we are matching each names(.PATTERN.TO.CV) to `pattern`, not the
    ## other way round. So we can provide a string describing the processing
    ## step and return the best matching CV terms.
    matches <- vapply(names(.PATTERN.TO.CV), FUN = function(z) {
        length(grep(pattern = z, x = pattern, ignore.case = TRUE)) > 0
    }, FUN.VALUE = logical(1))
    if (any(matches))
        res <- .PATTERN.TO.CV[matches]
    unname(res)
}
