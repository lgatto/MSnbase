#' @rdname Chromatograms-class
#'
#' @param object For \code{show}: a \code{Chromatograms} object.
setMethod("show", "Chromatograms", function(object) {
    nr <- nrow(object)
    nc <- ncol(object)
    cat(class(object), " with ",
        nr, ifelse(nr == 1, " row and ", " rows and "),
        nc, ifelse(nc == 1, " column\n", " columns\n"),
        sep = "")
    if (nr > 0 && nc > 0) {
        nms <- rownames(object)
        out <- apply(object, MARGIN = c(1, 2), function(z) {
            paste0("Chromatogram (", length(z[[1]]), ")")
        })
        if (!is.null(nms))
            rownames(out) <- nms
        print(out, quote = FALSE, right = TRUE)
    }
})

#' @rdname Chromatograms-class
#'
#' @description \code{Chromatograms} objects can, just like a \code{matrix},
#'     be subsetted using the \code{[} method. Single elements, rows or columns
#'     can be replaced using e.g. \code{x[1, 1] <- value} where \code{value}
#'     has to be a \code{Chromatogram} object or a \code{list} of such objects.
#'
#' @param x For \code{[}: the \code{Chromatograms} object to subset.
#' 
#' @param i For \code{[}: \code{numeric}, \code{logical} or \code{character}
#'     defining which row(s) to extract.
#'
#' @param j For \code{[}: \code{numeric}, \code{logical} or \code{character}
#'     defining which columns(s) to extract.
#'
#' @param drop For \code{[}: \code{logical(1)} whether to drop the
#'     dimensionality of the returned object (if possible).
#'
#' @return For \code{[}: the subset of the \code{Chromatograms} object. If a
#'     single element is extracted (e.g. if \code{i} and \code{j} are of length
#'     1) a \code{\link{Chromatogram}} object is returned. Otherwise a
#'     \code{list} of \code{\link{Chromatogram}} objects is returned, except
#'     if \code{drop = FALSE} in which case \emph{always} a \code{Chromatograms}
#'     object is returned.
setMethod("[", "Chromatograms",
          function(x, i, j, drop = TRUE) {
              if (missing(i) & missing(j))
                  return(x)
              if (missing(i))
                  i <- seq_len(nrow(x))
              if (missing(j))
                  j <- seq_len(ncol(x))
              if (is.logical(i))
                  i <- which(i)
              if (is.logical(j))
                  j <- which(j)
              x <- x@.Data[i = i, j = j, drop = drop]
              if (length(i) == 1 & length(j) == 1 & drop)
                  x <- x[[1]]
              if (!drop)
                  x <- as(x, "Chromatograms")
              if (validObject(x))
              x
          })

#' @rdname Chromatograms-class
#'
#' @param value For \code{[<-}: the replacement object(s). Can be a \code{list}
#'     of \code{\link{Chromatogram}} objects or, if length of \code{i} and
#'     \code{j} are 1, a single \code{\link{Chromatogram}} object.
setReplaceMethod("[", "Chromatograms",
                 function(x, i, j, value) {
                     if(missing(i) & missing(j))
                         return(x)
                     if (missing(i))
                         i <- seq_len(nrow(x))
                     if (missing(j))
                         j <- seq_len(ncol(x))
                     if (is.logical(i))
                         i <- which(i)
                     if (is.logical(j))
                         j <- which(j)
                     if (is.character(i))
                         stop("Sub-setting rows by name currently not supported")
                     if (is.character(j))
                         j <- match(j, colnames(x))
                     ## check i and j - currently we support only replacing
                     ## single elements, rows or columns.
                     if (length(i) > 1 & length(j) > 1)
                         stop("Only replacement of single elements, rows or ",
                              "columns supported")
                     ## Check value.
                     if (is(value, "Chromatogram"))
                         value <- list(value)
                     if (!all(sapply(value, function(z) is(z, "Chromatogram"))))
                         stop("All elements in 'value' have to be of type ",
                              "'Chromatogram'")
                     x@.Data[i, j] <- value
                     if (validObject(x))
                         x
                 })
