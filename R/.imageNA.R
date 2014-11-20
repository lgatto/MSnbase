##' @title NA heatmap visualisation for 2 groups
##' @param object An instance of class MSnSet
##' @param pcol Either the name of a phenoData variable to be used to
##' determine the group structure or a factor or any object that can
##' be coerced as a factor of length equal to nrow(object). The
##' resulting factor must have 2 levels. If missing (default)
##' \code{image(object)} is called.
##' @param Rowv Determines if and how the rows/features are
##' reordered. If missing (default), rows are reordered according to
##' \code{order((nNA1 + 1)^2/(nNA2 + 1))}, where NA1 and NA2 are the
##' number of missing values in each group. Use a vector of numerics
##' of feautre names to customise row order.
##' @param Colv A \code{logical} that determines if columns/samples
##' are reordered. Default is \code{TRUE}.
##' @param useGroupMean Replace individual feature intensities by the
##' group mean intensity. Default is FALSE.
##' @param ... Additional arguments passed to \code{\link{image}}.
##' @return Used for its side effect of plotting
##' @author Laurent Gatto
.imageNA2 <- function(object, pcol,
                      Rowv, Colv = TRUE,
                      useGroupMean = FALSE, ...) {   
    if (missing(pcol))
        image(object)
    if (is.character(pcol)) {
        if (pcol %in% varLabels(phenoData(object)))
            pcol <- factor(pData(object)[, pcol])
        else
            stop(pcol, " not found in varLabels(phenoData(object)): ",
                 paste(varLabels(phenoData(object)), collapse = " "))
    }
    pcol <- as.factor(pcol)
    stopifnot(length(levels(pcol)) == 2)
    g1 <- pcol == levels(pcol)[1]
    g2 <- pcol == levels(pcol)[2]        
    if (missing(Rowv)) {
        nNA1 <- apply(exprs(object), 1, function(x) sum(is.na(x[g1])))
        nNA2 <- apply(exprs(xx), 1, function(x) sum(is.na(x[g2])))
        Rowv <- order((nNA1 + 1)^2/(nNA2 + 1))
    }
    fn0 <- featureNames(object)
    object <- object[Rowv, ]
    if (identical(sort(fn0), featureNames(object)))
        warning("Feature names are different before and after reordering.")
    if (Colv) {
        ## reordering each protein values individually
        for (i in 1:nrow(object)) {
            k <- exprs(object)[i, ]
            k[rev(which(g1))] <- k[g1][order(k[g1])]
            k[g2] <- k[g2][order(k[g2])]
            exprs(object)[i, ] <- k   
        }
    }    
    if (useGroupMean) {
        for (i in 1:nrow(object)) {
            k1 <- exprs(object)[i, g1]        
            m1 <- mean(k1, na.rm=TRUE)    
            k1[!is.na(k1)] <- m1    
            k1[rev(which(g1))] <- k1[order(k1)]
            k2 <- exprs(object)[i, g2]       
            m1 <- mean(k2, na.rm=TRUE)    
            k2[!is.na(k2)] <- m1    
            k2 <- k2[order(k2)]
            exprs(object)[i, ] <- c(k1, k2)
        }        
    }
    image(object, ...)    
}
