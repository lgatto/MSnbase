## functions for Chromatograms class

.validChromatograms <- function(x) {
    msg <- character()
    ## All elements have to be of type Chromatogram
    if (length(x)) {
        res <- sapply(x, FUN = function(z) is(z, "Chromatogram"))
        if (!all(res))
            msg <- c(msg, paste0("All elements have to be of type ",
                                 "'Chromatogram'."))
        ## Shall we also ensure that fromFile in each column is the same?
    }
    if (length(msg))
        msg
    else TRUE
}

#' @description \code{Chromatograms}: create an instance of class
#'     \code{Chromatograms}.
#'
#' @param data A \code{list} of \code{\link{Chromatogram}} objects.
#'
#' @param ... Additional parameters to be passed to the \code{\link{matrix}}
#'     constructor, such as \code{nrow}, \code{ncol} and \code{byrow}.
#' 
#' @rdname Chromatograms-class
Chromatograms <- function(data, ...) {
    if (missing(data))
        return(new("Chromatograms"))
    res <- as(matrix(data, ...), "Chromatograms")
    if (validObject(res))
        res
}
