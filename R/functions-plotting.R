#' Takes the values for a single file.
#'
#' @param x `data.frame` with columns `"mz"`, `"rt"` and `"i"`.
#'
#' @param main `character(1)` with the title of the plot.
#'
#' @param col color for the circles.
#'
#' @param colramp color ramp to be used for the points' background.
#'
#' @param grid.color color to be used for the grid lines (or `NA` if they should
#'     not be plotted.
#'
#' @param pch The plotting character.
#'
#' @param layout `matrix` defining the layout of the plot, or `NULL` if
#'     `layout` was already called.
#'
#' @param ... additional parameters to be passed to the `plot` function.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @noRd
.plotXIC <- function(x, main = "", col = "grey", colramp = topo.colors,
                     grid.color = "lightgrey", pch = 21,
                     layout = matrix(1:2, ncol = 1), ...) {
    if (is.matrix(layout))
        layout(layout)
    if (nrow(x)) {
        ## Chromatogram.
        bpi <- unlist(lapply(split(x$i, x$rt), max, na.rm = TRUE))
        brks <- do.breaks(range(x$i), nint = 256)
        par(mar = c(0, 4, 2, 1))
        plot(as.numeric(names(bpi)), bpi, xaxt = "n", col = col, main = main,
             bg = level.colors(bpi, at = brks, col.regions = colramp),
             xlab = "", pch = pch, ylab = "", las = 2, ...)
        mtext(side = 4, line = 0, "Intensity", cex = par("cex.lab"))
        grid(col = grid.color)
        par(mar = c(3.5, 4, 0, 1))
        plot(x$rt, x$mz, main = "", pch = pch, col = col, xlab = "", ylab = "",
             yaxt = "n", bg = level.colors(x$i, at = brks,
                                           col.regions = colramp), ...)
        axis(side = 2, las = 2)
    } else {
        par(mar = c(0, 4, 2, 1))
        plot(3, 3, pch = NA, main = main, xlab = "", ylab = "", xaxt = "n",
             yaxt = "n", ...)
        mtext(side = 4, line = 0, "Intensity", cex = par("cex.lab"))
        grid(col = grid.color)
        par(mar = c(3.5, 4, 0, 1))
        plot(3, 3, pch = NA, main = "", xlab = "", ylab = "",
             xaxt = "n", yaxt = "n", ...)
    }
    grid(col = grid.color)
    mtext(side = 1, line = 2.5, "Retention time", cex = par("cex.lab"))
    mtext(side = 4, line = 0, "m/z", cex = par("cex.lab"))
}

#' Create a `matrix` to be used for the `layout` function to allow plotting of
#' vertically arranged *sub-plots* consisting of `sub_plot` plots.
#'
#' @param x `integer(1)` with the number of sub-plots.
#'
#' @param sub_plot `integer(1)` with the number of sub-plots per cell/plot.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' ## Assum we've got 5 *features* to plot and we want to have a two plots for
#' ## each feature arranged below each other.
#'
#' .vertical_sub_layout(5, sub_plot = 2)
.vertical_sub_layout <- function(x, sub_plot = 2) {
    sqrt_x <- sqrt(x)
    ncol <- ceiling(sqrt_x)
    nrow <- round(sqrt_x)
    rws <- split(1:(ncol * nrow * sub_plot), f = rep(1:nrow,
                                                     each = sub_plot * ncol))
    do.call(rbind, lapply(rws, matrix, ncol = ncol))
}
