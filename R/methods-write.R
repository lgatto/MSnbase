## write methods.

#' @rdname writeMSData
setMethod("write", signature(x = "ANY"), function(x, ...) {
    base::write(x, ...)
})

#' @title Write MS data to mzML or mzXML files
#'
#' @aliases write
#' 
#' @description The `write,MSnExp` and `write,OnDiskMSnExp` saves the content 
#'     of a [MSnExp] or [OnDiskMSnExp] object to MS file(s) in either
#'     *mzML* or *mzXML* format (one file per sample).
#'
#' @details The `write` method uses the *proteowizard* libraries through the
#'     `mzR` package to save the MS data. The data can be written to *mzML* or
#'     *mzXML* files with or without copying additional metadata information
#'     from the original files from which the data was read by the
#'     [readMSData()] function. This can be set using the `copy` parameter.
#'     Note that `copy = TRUE` requires the original files to be available.
#'
#'     The `write,ANY` method aims to restore the use of the [base::write()]
#'     function if `x` is neither a `MSnExp` nor an `OnDiskMSnExp` object.
#' 
#' @note General spectrum data such as total ion current, peak count, base peak
#'     m/z or base peak intensity are calculated from the actual spectrum data
#'     before writing the data to the files.
#'
#'     For MSn data, if the `OnDiskMSnExp` or `MSnExp` does not contain also the
#'     precursor scan of a MS level > 1 spectrum `precursorMZ`,
#'     `precursorCharge` and `precursorIntensity` of the spectrum is set to 0
#'     in the output file.
#' 
#' @param x `OnDiskMSnExp` or `MSnExp` object.
#'
#' @param files `character` with the file name(s). Its length has to match the
#'     number of samples/files of `x`.
#'
#' @param outformat `character(1)` defining the format of the output files.
#'     Default output format is `"mzml"`.
#'
#' @param verbose `logical(1)` if progress messages should be displayed.
#'
#' @param copy `logical(1)` if metadata (data processings, original file names
#'     etc) should be copied from the original files. See details for more
#'     information.
#' 
#' @param software_processing optionally provide specific data processing steps.
#'     See documentation of the `software_processing` parameter of
#'     [mzR::writeMSData()].
#'
#' @param ... For `write,ANY`: additional arguments for the [base::write()]
#'     function.
#' 
#' @author Johannes Rainer
#' 
#' @md
#'
#' @rdname writeMSData
setMethod("write", signature(x = "MSnExp"),
          function(x, files, outformat = c("mzml", "mzxml"),
                   verbose = isMSnbaseVerbose(), copy = TRUE,
                   software_processing = NULL) {
              ## Set copy to false if not all original files are available.
              if (copy & !all(file.exists(fileNames(x)))) {
                  warning("Setting 'copy = FALSE' because the original files ",
                          "can not be found.")
                  copy <- FALSE
              }
              ## Have to generate the header data.frame from the spectra.
              .writeMSData(x = x, files = files, outformat = outformat,
                           verbose = verbose, copy = copy, software_processing)
          })

#' @rdname writeMSData
setMethod("write", signature(x = "OnDiskMSnExp"),
          function(x, files, outformat = c("mzml", "mzxml"),
                   verbose = isMSnbaseVerbose(), copy = TRUE,
                   software_processing = NULL) {
              .writeMSData(x = x, files = files, outformat = outformat,
                           verbose = verbose, copy = copy, software_processing)
})

