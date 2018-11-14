plot_Spectrum1 <- function(spectrum,
                           centroided. = centroided(spectrum),
                           plot = TRUE) {
    if (is.na(centroided.)) {
        centroided. <- TRUE
        message("No centroided/profile information. Setting to centroided")
    }
    mtc <- mz(spectrum)
    i <- intensity(spectrum)
    dfr <- data.frame(i = i, mz = mtc)
    if (centroided.) {
        p <- ggplot(dfr, aes(x = mz, xend = mz, y = 0, yend = i)) +
            geom_segment()
    } else {
        p <- ggplot(dfr, aes(x = mtc, y = i)) + geom_line()
    }
    title <- ggtitle(paste("Retention time", rtime(spectrum)))
    p <- p + labs(x = "M/Z", y = "Intensity") + title
    if (plot)
        print(p + title)
    invisible(p + title)
}
