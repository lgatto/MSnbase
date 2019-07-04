## set plotNA method
.preparePlotNAData <- function(x) {
  ## pNA: percentrage of NAs allowed per feature
  nnacol <- colSums(is.na(x))
  nnarow <- rowSums(is.na(x))

  ocol <- order(nnacol)
  orow <- order(nnarow)

  nnacol <- nnacol[ocol]
  nnarow <- nnarow[orow]
  x <- x[orow, ocol]

  nc <- ncol(x)

  ## percentage of NA for each protein
  p <- 1 - nnarow / nc
  ## percentage of NA in data set
  d <- 1 - cumsum(nnarow) / (1:nrow(x) * nc)

  data.frame(x = seq_along(p), proteins = p, data = d)
}

plotNA_matrix <- function(X, pNA) {
  pNA <- pNA[1]
  dfr1 <- .preparePlotNAData(X)
  dfr2 <- data.frame(x = 1:nrow(dfr1),
                     variable = rep(c("proteins", "data"), each = nrow(dfr1)),
                     value = c(dfr1$proteins, dfr1$data))
  nkeep <- sum(dfr1$proteins >= (1 - pNA))
  kkeep <- dfr1$data[nkeep]
  x <- y <- z <- value <- variable <- NULL

  p <- ggplot() +
    geom_line(data = dfr2, aes(x = x, y = value, colour = variable)) +
      labs(x = "Protein index (ordered by data completeness)",
           y = "Data completeness") +
             theme(legend.position=c(0.23, 0.18),
                  legend.title = element_blank(),
                  legend.background = element_rect(size = 0)) +
                    scale_colour_hue(labels = c("Individual features", "Full dataset"),
                                     breaks = c("proteins", "data"))
  dfr0 <- data.frame(x = nrow(dfr1), y = min(dfr1$data))  
  p <- p +
    geom_point(data = dfr0, aes(x = x, y = y), alpha = 1/3) + 
      geom_text(data = dfr0, 
                aes(x = x, y = y, label = round(y, 2)),
                vjust = 1.5, size = 2.5)                 
  ## p <- p +
  ##   geom_text(data = dfr1, 
  ##             aes(x = length(proteins), y = min(data), label = round(min(data), 2)),
  ##             vjust = 1.5, size = 2.5) +
  ##               geom_point(data = dfr1, 
  ##                          aes(x = length(proteins), y = min(data)), alpha = 1/3)
  p <- p + 
    geom_text(data = data.frame(x = nkeep, y = kkeep),
              aes(x = x, y = y, label = round(y, 2)),
              hjust = -0.5, vjust = -0.5, size = 2.5) +              
                geom_point(data = data.frame(x = nkeep, y = kkeep),
                           aes(x = x, y = y), alpha = 1/3)
  
  p <- p + geom_text(data = data.frame(x = nkeep, y = (1 - pNA), z = nkeep),
                     aes(x = x, y = y, label = round(z, 2)),
                     size = 2.5, vjust = 2, hjust = 2) +
                       geom_point(data = data.frame(x = nkeep, y = (1 - pNA)),
                                  aes(x = x, y = y), alpha = 1/3)
  
  p <- p + annotate("text", label = nrow(X), x = 0, y = 1,
                    size = 2.5, vjust = -1, alpha = 1/3)

  print(p)
  invisible(p)
}



##' Produces a heatmap after reordring rows and columsn to highlight
##' missing value patterns.
##'
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
##' @param plot A \code{logical} specifying of an image should be
##' produced. Default is \code{TRUE}.
##' @param ... Additional arguments passed to \code{\link{image}}.
##' @return Used for its side effect of plotting. Invisibly returns
##' Rovw and Colv.
##' @author Laurent Gatto, Samuel Wieczorek and Thomas Burger
##' @examples
##' library("pRolocdata")
##' library("pRoloc")
##' data(dunkley2006)
##' pcol <- ifelse(dunkley2006$fraction <= 5, "A", "B")
##' nax <- makeNaData(dunkley2006, pNA = 0.10)
##' exprs(nax)[sample(nrow(nax), 30), pcol == "A"] <- NA
##' exprs(nax)[sample(nrow(nax), 50), pcol == "B"] <- NA
##' MSnbase:::imageNA2(nax, pcol)
##' MSnbase:::imageNA2(nax, pcol, useGroupMean = TRUE)
##' MSnbase:::imageNA2(nax, pcol, Colv = FALSE, useGroupMean = FALSE)
##' MSnbase:::imageNA2(nax, pcol, Colv = FALSE, useGroupMean = TRUE)
imageNA2 <- function(object, pcol,
                     Rowv, Colv = TRUE,
                     useGroupMean = FALSE,
                     plot = TRUE,
                     ...) {
    if (missing(pcol)) {
        if (plot) image2(object)
        return(invisible(NULL))
    }
    if (is.character(pcol) & length(pcol) == 1) {
        if (pcol %in% varLabels(phenoData(object))) {
            pcol <- factor(pData(object)[, pcol])
        }
        else
            stop(pcol, " not found in varLabels(phenoData(object)): ",
                 paste(varLabels(phenoData(object)), collapse = " "))
    }
    pcol <- as.factor(pcol)
    stopifnot(length(levels(pcol)) == 2)
    po <- order(pcol)
    pcol <- pcol[po]
    object <- object[, po]    
    g1 <- pcol == levels(pcol)[1]
    g2 <- pcol == levels(pcol)[2]        
    if (missing(Rowv)) {
        nNA1 <- apply(exprs(object), 1, function(x) sum(is.na(x[g1])))
        nNA2 <- apply(exprs(object), 1, function(x) sum(is.na(x[g2])))
        Rowv <- order((nNA1 + 1)^2/(nNA2 + 1))
    }
    fn0 <- featureNames(object)
    object <- object[Rowv, ]
    if (!identical(sort(fn0), sort(featureNames(object))))
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
            k2 <- exprs(object)[i, g2]       
            m1 <- mean(k2, na.rm=TRUE)    
            k2[!is.na(k2)] <- m1    
            exprs(object)[i, ] <- c(k1, k2)
        }        
    }
    if (plot) image2(object, ...)
    invisible(Rowv)
}


##' Visualise missing values as a heatmap and barplots along the
##' samples and features.
##'
##' @title Overview of missing value
##' @param object An object of class \code{MSnSet}.
##' @param verbose If verbose (default is \code{isMSnbaseVerbose()}), print a
##'     table of missing values.
##' @param reorderRows If reorderRows (default is \code{TRUE}) rows are ordered
##'     by number of NA.
##' @param reorderColumns If reorderColumns (default is \code{TRUE}) columns
##'     are ordered by number of NA.
##' @param ... Additional parameters passed to \code{image2}.
##' @return Used for its side effect. Invisibly returns \code{NULL}
##' @author Laurent Gatto
##' @examples
##' data(naset)
##' naplot(naset)
naplot <- function(object, verbose=isMSnbaseVerbose(),
                   reorderRows=TRUE, reorderColumns=TRUE, ...) {
    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    zones <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(zones, widths=c(4/5, 1/5), heights=c(1/5, 4/5))
    mNA <- is.na(exprs(object))
    features.na <- rowSums(mNA)
    samples.na <- colSums(mNA)
    if (reorderRows) {
      xo <- order(features.na)
    } else {
      xo <- 1L:nrow(object)
    }
    if (reorderColumns) {
      yo <- order(samples.na)
    } else {
      yo <- 1L:ncol(object)
    }
    par(mar=c(3,3,1,1))
    image2(object[xo, yo], ...)
    par(mar=c(0,3,1,1))
    barplot(samples.na[yo], space=0, xaxt="n", xaxs="i")
    par(mar=c(3,0,1,1))
    barplot(features.na[xo], space=0, horiz=TRUE, yaxt="n", yaxs="i")
    if (verbose) {
        print(table(features.na))
        print(table(samples.na))
    }
    invisible(NULL)
}
