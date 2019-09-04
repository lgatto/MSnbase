##' plot spectrum1 vs spectrum2
##' @param spectra list, 2 MSnbase::Spectrum2 objects
##' @param tolerance double, allowed deviation to be considered as equal peaks
##' @param relative relative (or absolute) deviation
##' @param ... additional paramters passed to \code{.plotSpectrumVsSpectrum}.
##' @noRd
plotSpectrumVsSpectrum <- function(spectra, tolerance=25e-6,
                                   relative=TRUE,
                                   ...) {
  common <- lapply(list(c(1, 2), c(2, 1)), function(x) {
    commonPeaks(spectra[[x[1]]], spectra[[x[2]]],
                method="highest", tolerance=tolerance, relative=relative)
  })
  .plotSpectrumVsSpectrum(spectra, common=common, ...)
}

#' plot spectrum1 vs spectrum2
##' @param spectra list, 2 MSnbase::Spectrum2 objects
##' @param sequences character vector (length==2) containing the peptide
##' sequences for both spectra
##' @param common list (length==2), containing logical vector for common peaks
##' @param norm normalize?
##' @param xlim limits for x-axis
##' @param ylim limits for y-axis
##' @param legend.cex cex for legend
##' @param peaks.pch pch for marking peaks
##' @param ... additional parameters passed to \code{.plotSingleSpectrum}.
##' @param fragments.cex cex for fragments
##' @noRd
.plotSpectrumVsSpectrum <- function(spectra,
                                    sequences,
                                    common,
                                    norm = TRUE,
                                    xlim, ylim,
                                    legend.cex = 1,
                                    peaks.pch = 19, ...) {
  if (norm) 
    spectra <- lapply(spectra, normalize)

  if (missing(xlim)) {
    mass <- unlist(lapply(spectra, mz))
    xlim <- c(min(c(mass, 0)), max(c(mass, 0)))
  }

  if (missing(ylim)) {
    inten <- unlist(lapply(spectra, intensity))
    maxInten <- max(c(inten, 0))
    ylim <- c(-maxInten, maxInten)
  }

  if (missing(common)) 
    common <- lapply(spectra, function(x)logical(peaksCount(x)))

  if (missing(sequences)) 
    sequences <- character(2)

  orientation <- c(1, -1)
  add <- c(FALSE, TRUE)
  legend.pos <- c("topleft", "bottomleft")
  ## colors: ColorBrewer RdYlBu c(9, 11, 3, 1)
  cols <- c("#74ADD1", "#313695", "#F46D43", "#A50026")
  pch <- c(NA, peaks.pch)

  for (i in seq(along = spectra)) {
    .plotSingleSpectrum(spectra[[i]], sequence=sequences[[i]],
                        orientation=orientation[i], add=add[i],
                        xlim=xlim, ylim=ylim,
                        col=cols[(i-1)*2+common[[i]]+1],
                        pch=pch[common[[i]]+1], ...)

    if (msLevel(spectra[[i]]) == 1) {
        label <- paste0("Retention time: ", formatRt(rtime(spectra[[i]])),
                        ", # common: ", sum(common[[i]]))
        
    } else {
        label <- paste0("prec scan: ", precScanNum(spectra[[i]]))
        if (peaksCount(spectra[[i]])) {
            label <- paste0(label, ", prec mass: ", round(precursorMz(spectra[[i]]), 3),
                            ", prec z: ", precursorCharge(spectra[[i]]),
                            ", # common: ", sum(common[[i]]))
            if (nchar(sequences[[i]])) {
                label <- paste0(label, ", seq: ", sequences[[i]])
            }
        }
    }
    legend(legend.pos[i], legend = label, bty="n", cex = legend.cex)
  }
}

#' plot spectrum
##' @param object MSnbase::Spectrum
##' @param sequence character, peptide sequence
##' @param orientation c(1,-1) up/down
##' @param add if add TRUE "plot" is not called and the spectrum is drawn in the
##' existing plot (see ?hist for a similar argument)
##' @param col col would be recycled to the length of peaksCount(object)
##' @param pch pch would be recycled to the length of peaksCount(object)
##' @param xlab label for the x-axis
##' @param ylab label for the y-axis
##' @param xlim limits for the x-axis
##' @param ylim limits for the y-axis
##' @param tolerance double, allowed deviation
##' @param relative relative (or absolute) deviation
##' @param type fragment types, could be c("a", "b", "c", "x", "y", "z")
##' @param modifications a named (amino acid one-letter-code; upper case) vector
##' @param neutralLoss list, currently water and ammonia loss are supported
##' @param z fragment charge
##' @param fragments a data.frame produced by calculatedFragments_Spectrum2
##' @param fragments.cex cex for the fragment letters
##' @param peaks.lwd lwd for the peaks
##' @param peaks.cex cex for the points of on the top of the peaks
##' @param ... further arguments passed to plot.default
##' @noRd
.plotSingleSpectrum <- function(object, sequence,
                                orientation = 1, add = FALSE,
                                col = "#74ADD1", pch = NA,
                                xlab = "m/z", ylab = "intensity",
                                xlim, ylim = c(0, 1),
                                tolerance = 0.1, relative = FALSE,
                                type = c("b", "y"),
                                modifications = c(C = 57.02146),
                                neutralLoss = defaultNeutralLoss(),
                                z = 1,
                                fragments = calculateFragments_Spectrum2(object,
                                  sequence = sequence, tolerance = tolerance,
                                  relative = relative, type = type, z = z,
                                  modifications = modifications,
                                  neutralLoss = neutralLoss,
                                  verbose = isMSnbaseVerbose()),
                               fragments.cex = 0.75, peaks.lwd = 1, peaks.cex = 0.5, ...) {
  if (peaksCount(object) > 0 && !centroided(object)) {
    message("Your spectrum is not centroided.")
  }

  if (missing(xlim)) {
    xlim <- c(min(c(mz(object), 0)), max(c(mz(object), 0)))
  }

  if (missing(ylim)) {
    ylim <- c(0, max(c(intensity(object), 0)))
  }

  if (!add) {
    plot(NA, type = "h", col = 1,
         xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, ...)
    abline(h = 0, col = "#808080")
  }

  lines(mz(object), orientation*intensity(object),
        type="h", col=col, lwd=peaks.lwd)
  points(mz(object), orientation*intensity(object),
         col=col, pch=pch, cex=peaks.cex)

  if (nrow(fragments)) {
    text(fragments$mz,
         orientation*intensity(object)[match(fragments$mz, mz(object))],
         fragments$ion, pos=2+orientation, offset=0.25,
         cex=fragments.cex, col="#808080")
  }
}

