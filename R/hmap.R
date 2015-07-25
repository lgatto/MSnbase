## hmap <- function(x, fcol, pcol, ...) {
##     hap <- new("HeatmapAnnotation")
##     if (missing(pcol) & ncol(pData(x)) > 0) {
##         hap <- HeatmapAnnotation(pData(x))
##     } else {
##         hap <- HeatmapAnnotation(pData(x)[, pcol, drop=FALSE])
##     }
##     hm <- Heatmap(exprs(x), top_annotation = hap, ...)
##     if (!missing(fcol)) {
##         haf <- rowAnnotation(fData(x)[, fcol, drop = FALSE])
##         hm <- hm + haf
##     }
##     print(hm)
##     invisible(hm)
## }
