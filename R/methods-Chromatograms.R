#' @rdname Chromatograms-class
#'
#' @param object a \code{Chromatograms} object.
setMethod("show", "Chromatograms", function(object) {
    nr <- nrow(object)
    nc <- ncol(object)
    cat(class(object), " with ",
        nr, ifelse(nr == 1, " row and ", " rows and "),
        nc, ifelse(nc == 1, " column\n", " columns\n"),
        sep = "")
    sumFun <- function(z) {
        paste0("length: ", length(z[[1]]))
    }
    if (nr > 0 && nc > 0) {
        if (nr <= 4) {
            out <- apply(object, MARGIN = c(1, 2), sumFun)
            rownames(out) <- paste0("[", 1:nrow(out), ",]")
        }
        else {
            out <- rbind(
                apply(object[c(1, 2), , drop = FALSE], MARGIN = c(1, 2), sumFun),
                rep(" ... ", ncol(object)),
                apply(object[nrow(object) - c(1, 0), , drop = FALSE],
                      MARGIN = c(1, 2), sumFun)
            )
            rownames(out) <- c("[1,]", "[2,]", "...",
                               paste0("[", c(nrow(object) - c(1, 0)), ",]"))
        }
        rn <- rownames(out)
        out <- rbind(rep("<Chromatogram>", ncol(out)), out)
        rownames(out) <- c("", rn)        
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
#' @note Subsetting with \code{[} will always return a \code{Chromatograms}
#'     object (with the exception of extracting a single element)
#'     unless \code{drop = TRUE} is specified. This is different from the
#'     default subsetting behaviour of \code{matrix}-like objects.
#' 
#' @param x For all methods: a \code{Chromatograms} object.
#' 
#' @param i For \code{[}: \code{numeric}, \code{logical} or \code{character}
#'     defining which row(s) to extract.
#'
#' @param j For \code{[}: \code{numeric}, \code{logical} or \code{character}
#'     defining which columns(s) to extract.
#'
#' @param drop For \code{[}: \code{logical(1)} whether to drop the
#'     dimensionality of the returned object (if possible). The default is
#'     \code{drop = FALSE}, i.e. each subsetting returns a \code{Chromatograms}
#'     object (or a \code{Chromatogram} object if a single element is
#'     extracted).
#'
#' @return For \code{[}: the subset of the \code{Chromatograms} object. If a
#'     single element is extracted (e.g. if \code{i} and \code{j} are of length
#'     1) a \code{\link{Chromatogram}} object is returned. Otherwise (if
#'     \code{drop = FALSE}, the default, is specified) a \code{Chromatograms}
#'     object is returned. If \code{drop = TRUE} is specified, the method
#'     returns a \code{list} of \code{Chromatogram} objects.
#'
#'     For \code{phenoData}: an \code{NAnnotatedDataFrame} representing the
#'     pheno data of the object.
#'
#'     For \code{pData}: a \code{data.frame} representing the pheno data of
#'     the object.
#'
#'     For \code{$}: the value of the corresponding column in the pheno data
#'     table of the object.
setMethod("[", "Chromatograms",
          function(x, i, j, drop = FALSE) {
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
              ## Return a single element as a Chromatogram
              if (length(i) == 1 & length(j) == 1)
                  return(x@.Data[i, j, drop = TRUE][[1]])
              pd <- x@phenoData
              ## Multiple elements, return type depends on drop.
              x <- x@.Data[i = i, j = j, drop = drop]
              if (!drop) {
                  x <- as(x, "Chromatograms")
                  pd <- pd[j, ]
                  ## Drop levels
                  pData(pd) <- droplevels(pData(pd))
                  ## set row names
                  rownames(pd) <- colnames(x)
                  x@phenoData <- pd
              }
              if (validObject(x))
                  x
          })

#' @rdname Chromatograms-class
#'
#' @param value For \code{[<-}: the replacement object(s). Can be a \code{list}
#'     of \code{\link{Chromatogram}} objects or, if length of \code{i} and
#'     \code{j} are 1, a single \code{\link{Chromatogram}} object.
#' 
#'     For \code{pData<-}: a \code{data.frame} with the number of rows matching
#'     the number of columns of \code{object}.
#'
#'     For \code{colnames}: a \code{character} with the new column names.
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
                     res <- vapply(value, FUN = is, FUN.VALUE = logical(1),
                                   "Chromatogram")
                     if(!all(res))
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
                      .plotChromatogramList(x[i, , drop = TRUE], col = col,
                                            lty = lty, type = type, xlab = xlab,
                                            ylab = ylab, main = main, ...)
                  else
                      plot(x[i, 1], col = col, lty = lty, type = type,
                           xlab = xlab, ylab = ylab, main = main, ...)
              }
          })

#' @rdname Chromatograms-class
#'
#' @description \code{phenoData}: accesses the phenotypical desccription of the
#'     samples. Returns an \code{NAnnotatedDataFrame} object.
setMethod("phenoData", "Chromatograms", function(object) object@phenoData)

#' @rdname Chromatograms-class
#'
#' @description \code{pData}: accesses the phenotypical description of the
#'     samples. Returns a \code{data.frame}.
setMethod("pData", "Chromatograms", function(object) pData(phenoData(object)))

#' @rdname Chromatograms-class
#'
#' @description \code{pData<-}: replace the phenotype data.
setReplaceMethod("pData", "Chromatograms", function(object, value) {
    pData(object@phenoData) <- value
    if (any(rownames(value) != as.character(1:nrow(value))))
        colnames(object@.Data) <- rownames(value)
    if (validObject(object))
        object
})

#' @rdname Chromatograms-class
#'
#' @description \code{$} and \code{$<-}: get or replace individual columns of
#'     the object's pheno data.
#'
#' @param name For \code{$}, the name of the pheno data column.
setMethod("$", "Chromatograms", function(x, name) {
    eval(substitute(pData(x)$NAME_ARG, list(NAME_ARG = name)))
})
#' @rdname Chromatograms-class
setReplaceMethod("$", "Chromatograms", function(x, name, value) {
    pData(x)[[name]] <- value
    x
})

#' @rdname Chromatograms-class
#'
#' @description \code{colnames<-}: replace or set the column names of the
#'     \code{Chromatograms} object. Does also set the \code{rownames} of the
#'     \code{phenoData}.
setReplaceMethod("colnames", "Chromatograms", function(x, value) {
    rownames(pData(x)) <- value
    if (validObject(x))
        x
})
