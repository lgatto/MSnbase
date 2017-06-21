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

#' @rdname Chromatograms-class
#'
#' @description \code{plot}: plots a \code{Chromatograms} object. For each row
#'     in the object one plot is created, i.e. all \code{\link{Chromatogram}}
#'     objects in the same row are added to the same plot.
#'
#' @details \code{plot}: if \code{nrow(x) > 1} the plot area is split into
#'     \code{nrow(x)} sub-plots and the chromatograms of one row are plotted in
#'     each.
#'
#' @param col For \code{plot}: the color to be used for plotting. Either a
#'     vector of length 1 or equal to \code{ncol(x)}.
#'
#' @param lty For \code{plot}: the line type (see \code{\link[graphics]{plot}}
#'     for more details. Can be either a vector of length 1 or of length equal
#'     to \code{ncol(x)}.
#'
#' @param type For \code{plot}: the type of plot (see
#'     \code{\link[graphics]{plot}} for more details. Can be either a vector
#'     of length 1 or of length equal to \code{ncol(x)}.
#' 
#' @inheritParams Chromatogram-class
#'
#' @examples
#' 
#' ## Create some random Chromatogram objects
#' ints <- abs(rnorm(123, mean = 200, sd = 32))
#' ch1 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 231)
#' ints <- abs(rnorm(122, mean = 250, sd = 43))
#' ch2 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 231)
#' ints <- abs(rnorm(125, mean = 590, sd = 120))
#' ch3 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 542)
#' ints <- abs(rnorm(124, mean = 1200, sd = 509))
#' ch4 <- Chromatogram(rtime = seq_along(ints), intensity = ints, mz = 542)
#'
#' ## Combine into a 2x2 Chromatograms object
#' chrs <- Chromatograms(list(ch1, ch2, ch3, ch4), byrow = TRUE, ncol = 2)
#'
#' ## Plot the second row
#' plot(chrs[2, , drop = FALSE])
#'
#' ## Plot all chromatograms
#' plot(chrs, col = c("#ff000080", "#00ff0080"))
setMethod("plot", signature = signature("Chromatograms"),
          function(x, col = "#00000060", lty = 1, type = "l",
                   xlab = "retention time", ylab = "intensity",
                   main = NULL, ...){
              nr <- nrow(x)
              nc <- ncol(x)
              ## Initialize plot if we've got more than one row.
              if (nr > 1) {
                  par(mfrow = c(round(sqrt(nr)), ceiling(sqrt(nr))))
              }
              for (i in seq_len(nr)) {
                  if (nc > 1)
                      .plotChromatogramList(x[i, ], col = col, lty = lty,
                                            type = type, xlab = xlab,
                                            ylab = ylab, main = main, ...)
                  else
                      plot(x[i, 1], col = col, lty = lty, type = type,
                           xlab = xlab, ylab = ylab, main = main, ...)
              }
          })
