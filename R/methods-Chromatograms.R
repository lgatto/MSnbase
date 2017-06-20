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
