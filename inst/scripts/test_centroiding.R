## Testing centroiding with peak mz refining. The relevating code is in
## function-Spectrum.R (.refinePeakMz). All parameters to this function are
## passed via the ... parameters of the pickPeaks method.
library(MSnbase)
library(xcms)
library(RColorBrewer)

#' Create a plot that combines a XIC and a mz/rt 2D plot.
#'
#' @param x `data.frame` such as returned by the `extractMsData` function. Only
#'     a single `data.frame` has to be submitted.
#'
#' @param main `character(1)` specifying the title.
#'
#' @param mfrow `numeric(2)` defining the plot layout. This will be passed
#'     directly to `par(mfrow = mfrow)`. See `par` for more information. Setting
#'     `mfrow = NULL` avoids calling `par(mfrow = mfrow)` hence allowing to
#'     pre-define the plot layout.
#'
#' @author Johannes Rainer
#' 
#' @md
#' 
#' @noRd
plot_msdata <- function(x, main = "", cex = 1, mfrow = c(2, 1)) {    
    cols <- level.colors(x$i, at = do.breaks(range(x$i), nint = 256),
                         col.regions =
                             colorRampPalette(rev(brewer.pal(9, "YlGnBu"))))
    if (length(mfrow) == 2)
        par(mfrow = mfrow)
    par(mar = c(0, 4, 2, 1))
    x_split <- split(x, f = x$rt)
    ints <- lapply(x_split, function(z) sum(z$i))
    plot(as.numeric(names(ints)), ints, main = main, xlab = "", xaxt = "n",
         ylab = "", las = 2, pch = 21, bg = cols, col = "grey", cex = cex)
    mtext(side = 4, line = 0, "intensity", cex = par("cex.lab"))
    grid(col = "#00000020")
    par(mar = c(3.5, 4, 0, 1))
    plot(x$rt, x$mz, main = "", pch = 21, bg = cols, col = "grey",
         xlab = "", ylab = "", yaxt = "n", cex = cex)
    axis(side = 2, las = 2)
    grid(col = "#00000020")
    mtext(side = 1, line = 2.5, "retention time", cex = par("cex.lab"))
    mtext(side = 4, line = 0, "mz", cex = par("cex.lab"))
}

## The data files
mz_cent <- "/Users/jo/data/2016/2016_06/130616_POOL_IntraP_PC_POS_6.mzML"
mz_prof <- "/Users/jo/data/2016/2016_06_profile/130616_POOL_IntraP_PC_POS_6.mzML"
data_cent <- readMSData(mz_cent, mode = "onDisk")
data_prof <- readMSData(mz_prof, mode = "onDisk")

## 3-Methylhistidine
mz_exp <- 170.0918712
mzr <- c(170.0895, 170.0945)
rtr <- c(180, 194)

layout(matrix(1:6, nrow = 2))
## Read centroided data:
ms_cent <- extractMsData(data_cent, rt = rtr, mz = mzr)
plot_msdata(ms_cent[[1]], main = "proteowizard", mfrow = NULL)
abline(h = mz_exp, col = "red", lty = 3)

## Read profile data:
ms_prof <- extractMsData(data_prof, rt = rtr, mz = mzr)
plot_msdata(ms_prof[[1]], main = "profile", mfrow = NULL)
abline(h = mz_exp, col = "red", lty = 3)

## Standard peak picking
data_cent_R <- pickPeaks(data_prof)
ms_cent_R <- extractMsData(data_cent_R, rt = rtr, mz = mzr)
plot_msdata(ms_cent_R[[1]], main = "plain picking", mfrow = NULL)
abline(h = mz_exp, col = "red", lty = 3)

## Refine the peak based on 2 neighbors
data_cent_n2 <- pickPeaks(data_prof, refineMz = "neighbors", n = 2)
ms_cent_n2 <- extractMsData(data_cent_n2, rt = rtr, mz = mzr)
plot_msdata(ms_cent_n2[[1]], main = "2 neighbors", mfrow = NULL)
abline(h = mz_exp, col = "red", lty = 3)
