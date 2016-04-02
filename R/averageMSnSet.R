##' Given a list of \code{MSnSet} instances, typically representing
##' replicated experiments, the function returns an average
##' \code{MSnSet}.
##'
##' This function is aimed at facilitating the visualisation of
##' replicated experiments and should not be used as a replacement for
##' a statistical analysis.
##' 
##' The samples of the instances to be averaged must be identical but
##' can be in a different order (they will be reordered by
##' default). The features names of the result will correspond to the
##' union of the feature names of the input \code{MSnSet}
##' instances. Each average value will be computed by the \code{avg}
##' function and the dispersion of the replicated measurements will be
##' estimated by the \code{disp} function. These dispersions will be
##' stored as a \code{data.frame} in the feature metadata that can be
##' accessed with \code{fData(.)$disp}. Similarly, the number of
##' missing values that were present when average (and dispersion)
##' were computed are available in \code{fData(.)$disp}.
##'
##' Currently, the feature metadata of the returned object corresponds
##' the the feature metadata of the first object in the list
##' (augmented with the missing value and dispersion values); the
##' metadata of the features that were missing in this first input are
##' missing (i.e. populated with \code{NA}s). This may change in the
##' future.
##'
##' @title Generate an average \code{MSnSet}
##' @param x A \code{list} of valid \code{MSnSet} instances to be averaged.
##' @param avg The averaging function. Default is the mean after
##' removing missing values, as computed by \code{function(x)
##' mean(x, na.rm = TRUE)}.
##' @param disp The disperion function. Default is an non-parametric
##' coefficient of variation that replaces the standard deviation by
##' the median absolute deviation as computed by
##' \code{mad(x)/abs(mean(x))}. See \code{\link{npcv}} for
##' details. Note that the \code{mad} of a single value is 0 (as
##' opposed to \code{NA} for the standard deviation, see example
##' below).
##' @return A new average \code{MSnSet}.
##' @author Laurent Gatto
##' @seealso \code{\link{compfnames}} to compare MSnSet feature names.
##' @examples
##' library("pRolocdata")
##' ## 3 replicates from Tan et al. 2009
##' data(tan2009r1)
##' data(tan2009r2)
##' data(tan2009r3)
##' x <- MSnSetList(list(tan2009r1, tan2009r2, tan2009r3))
##' avg <- averageMSnSet(x)
##' dim(avg)
##' head(exprs(avg))
##' head(fData(avg)$nNA)
##' head(fData(avg)$disp)
##' ## using the standard deviation as measure of dispersion
##' avg2 <-averageMSnSet(x, disp = sd)
##' head(fData(avg2)$disp)
##' ## keep only complete observations, i.e proteins 
##' ## that had 0 missing values for all samples
##' sel <- apply(fData(avg)$nNA, 1 , function(x) all(x == 0))
##' avg <- avg[sel, ]
##' disp <- rowMax(fData(avg)$disp)
##' library("pRoloc")
##' setStockcol(paste0(getStockcol(), "AA"))
##' plot2D(avg, cex = 7.7 * disp)
##' title(main = paste("Dispersion: non-parametric CV",
##'                    paste(round(range(disp), 3), collapse = " - ")))
averageMSnSet <- function(x,
                          avg = function(x) mean(x, na.rm = TRUE),
                          disp = npcv) {
    if (!inherits(x, "MSnSetList"))
        stop("'x' must be an 'MSnSetList'")

    x <- msnsets(x)

    l <- length(x)
    if (l < 2) return(x[[1]])    
    if (is.null(names(x))) names(x) <- paste0("x", 1:l)

    ## reordering cols, just in case
    x <- lapply(x, function(.x) .x[, order(sampleNames(.x))])

    ## checking sample names
    sn <- lapply(x, sampleNames)
    if (!all(sapply(sn[-1], identical, sn[[1]])))
        stop("Sample names do not match.")

    m <- length(sn[[1]])
    fn <- unique(unlist(lapply(x, featureNames)))
    n <- length(fn)

    e <- array(NA, c(n, m, l))
    dimnames(e) <- list(fn, sn[[1]], names(x))

    for (i in 1:l) {
        k <- featureNames(x[[i]])
        e[k, , i] <- exprs(x[[i]])[k, ]
    }

    e0 <- apply(e, c(1, 2), avg)
    fd <- fData(x[[1]])
    fd0 <- data.frame(matrix(NA, ncol = ncol(fd), nrow = n - nrow(fd)))
    rownames(fd0) <- setdiff(fn, rownames(fd))
    names(fd0) <- names(fd)

    stopifnot(identical(rownames(e0)[1:nrow(fd)],
                        rownames(fd)))
    fd0 <- rbind(fd, fd0)
    stopifnot(all.equal(rownames(fd0), rownames(e0)))
    fd0$nNA <- apply(e, c(1, 2), function(i) sum(is.na(i)))
    fd0$disp <- apply(e, c(1, 2), disp)

    ans <- new("MSnSet",
               exprs = e0,
               featureData = new("AnnotatedDataFrame", data = fd0),
               phenoData = phenoData(x[[1]]))
    ans <- logging(ans, paste("Average of", l))
    if (validObject(ans)) ans
}
