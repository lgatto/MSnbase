\name{pickPeaks-methods}
\docType{methods}
\alias{pickPeaks-methods}
% \alias{pickPeaks,MSnExp-method}
% \alias{pickPeaks,Spectrum-method}
\alias{pickPeaks}


\title{ Peak Detection for 'MSnExp' or 'Spectrum' instances }

\description{
  This method performs a peak picking on individual spectra
  (\code{Spectrum} instances) or whole experiments (\code{MSnExp} instances) to
  create centroided spectra.
  For noisy spectra there are currently two different noise estimators, the
  Median Absolute Deviation (\code{method = "MAD"}) and
  Friedman's Super Smoother (\code{method = "SuperSmoother"}),
  as implemented in the \code{MALDIquant::detectPeaks} and
  \code{MALDIquant::estimateNoise} functions respectively.
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "MSnExp", halfWindowSize = "integer",
        method = "character", SNR = "numeric", verbose = "logical", ...)}}{
      Performs the peak picking for all spectra in an \code{MSnExp} instance.
      \code{method} could be \code{"MAD"} or \code{"SuperSmoother"}.
      \code{halfWindowSize} controls the window size of the peak picking
      algorithm. The resulting window size is \code{2 * halfWindowSize + 1}.
      The size should be nearly (or slightly larger) the \emph{FWHM}
      (full width at half maximum).
      A local maximum is considered as peak if its intensity is \code{SNR}
      times larger than the estimated noise.
      The arguments \code{\dots} are passed to the noise estimator functions.
      Currenlty only the \code{method = "SuperSmoother"} accepts additional
      arguments, e.g. \code{span}. Please see \code{\link[stats]{supsmu}} for
      details.
      This method displays a progress bar if \code{verbose = TRUE}.
      Returns an \code{MSnExp} instance with centroided spectra.
    }
    \item{\code{signature(x = "Spectrum", method = "character",
        halfWindowSize = "integer", ...)}}{
      Performs the peak picking for the spectrum (\code{Spectrum} instance).
      This method is the same as above but returns a centroided \code{Spectrum}
      instead of an \code{MSnExp} object. It has no \code{verbose} argument.
      Please read the details for the above \code{MSnExp} method.
    }
  }
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  \code{\link{clean}}, \code{\link{removePeaks}} \code{\link{smooth}},
  \code{\link{estimateNoise}} and \code{\link{trimMz}} for other spectra
  processing methods.
}

\references{
S. Gibb and K. Strimmer. 2012.
MALDIquant: a versatile R package for the analysis of mass spectrometry data.
Bioinformatics 28: 2270-2271.
\url{http://strimmerlab.org/software/maldiquant/}
}

\examples{
sp1 <- new("Spectrum1",
           intensity = c(1:6, 5:1),
           mz = 1:11,
           centroided = FALSE)
sp2 <- pickPeaks(sp1)
intensity(sp2)

data(itraqdata)
itraqdata2 <- pickPeaks(itraqdata)
processingData(itraqdata2)
}

\keyword{methods}
