\name{plot.Spectrum.Spectrum-methods}
\alias{plot.Spectrum.Spectrum}
\alias{plot,Spectrum,Spectrum-method}
\docType{methods}

\title{ Plotting a 'Spectrum' vs another 'Spectrum' object. }

\description{
  These method plot mass spectra MZ values against the intensities as
  line plots. The first spectrum is plotted in the upper panel and the other in
  upside down in the lower panel. Common peaks are drawn in a slightly darker
  colour. If a peptide sequence is provided it automatically calculates
  and labels the fragments.
}

\arguments{
  \item{x}{ Object of class \code{"\linkS4class{Spectrum}"} .}
  \item{y}{ Object of class \code{"\linkS4class{Spectrum}"} .}
  \item{\dots}{Further arguments passed to internal functions.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "Spectrum", y = "Spectrum", \dots)}}{
      Plots two spectra against each other. Common peaks are drawn in a slightly
      darker colour.
      The \code{\dots} arguments are passed to the internal functions.
      Currently \code{tolerance}, \code{relative}, \code{sequences} and most of
      the \code{\link{plot.default}} arguments (like \code{xlim}, \code{ylim},
      \code{main}, \code{xlab}, \code{ylab}, \dots) are supported.
      You could change the \code{tolerance} (default \code{25e-6}) and
      decide whether this tolerance should be applied relative
      (default \code{relative = TRUE}) or absolute (\code{relative = FALSE})
      to find and colour common peaks.
      Use a \code{character} vector of length 2 to provide \code{sequences}
      which would be used to calculate and draw the corresponding fragments.
      If \code{sequences} are given the
      \code{type} argument (default: \code{type=c("b", "y")} specify the
      fragment types which should calculated. Also it is possible to allow some
      \code{modifications}. Therefore you have to apply a named \code{character}
      vector for \code{modifications} where the name corresponds to the
      one-letter-code of the modified amino acid
      (default: Carbamidomethyl \code{modifications=c(C=57.02146)}). Additional
      you can specifiy the type of \code{neutralLoss} (default:
      \code{PSMatch::defaultNeutralLoss()}).
      See \code{\link{calculateFragments}} for details.

      There are a lot of graphical arguments available to control the
      representation of the peaks and fragments. Use \code{peaks.pch} to set the
      character on top of the peaks (default: \code{peaks.pch=19}). In a similar
      way you can set the line width \code{peaks.lwd=1} and the magnification
      \code{peaks.cex=0.5} of the peaks. The size of the fragment/legend labels could
      be set using \code{fragments.cex=0.75} or \code{legend.cex} respectively.
      See \code{\link{par}} for details about graphical parameters in general.
    }
  }
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  More spectrum plotting available in \code{\link{plot.Spectrum}}.

  More details about fragment calculation: \code{\link{calculateFragments}}.
}

\examples{
## find path to a mzXML file
file <- dir(system.file(package = "MSnbase", dir = "extdata"),
            full.name = TRUE, pattern = "mzXML$")

## create basic MSnExp
msexp <- readMSData(file, centroided.=FALSE)

## centroid them
msexp <- pickPeaks(msexp)

## plot the first against the second spectrum
plot(msexp[[1]], msexp[[2]])

## add sequence information
plot(msexp[[1]], msexp[[2]], sequences=c("VESITARHGEVLQLRPK",
                                         "IDGQWVTHQWLKK"))


itraqdata2 <- pickPeaks(itraqdata)
(k <- which(fData(itraqdata2)[, "PeptideSequence"] == "TAGIQIVADDLTVTNPK"))
mzk <- precursorMz(itraqdata2)[k]
zk <- precursorCharge(itraqdata2)[k]
mzk * zk
plot(itraqdata2[[k[1]]], itraqdata2[[k[2]]])
}

\keyword{methods}
