#' @noRd
.valid.BackendMemory.spectra.names <- function(x) {
    n <- length(x)
    nms <- names(x)

    if (n && any(is.null(nms)))
        return("Spectra names should not be NULL.")
    if (n && anyNA(nms))
        return("Spectra names should not contain NA.")
    if (n && !all(nchar(nms)))
        return("Spectra names should not be missing.")
    if (n && anyDuplicated(nms))
        return("Duplicated spectra names found.")
    NULL
}
