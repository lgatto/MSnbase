\name{estimateNoise-methods}
\docType{methods}
\alias{estimateNoise-methods}
\alias{estimateNoise}


\title{ Noise Estimation for 'Spectrum' instances }

\description{
  This method performs a noise estimation on individual spectra
  (\code{Spectrum} instances).
  There are currently two different noise estimators, the
  Median Absolute Deviation (\code{method = "MAD"}) and
  Friedman's Super Smoother (\code{method = "SuperSmoother"}),
  as implemented in the \code{MALDIquant::detectPeaks} and
  \code{MALDIquant::estimateNoise} functions respectively.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "Spectrum", method = "character", ...)}}{
      Estiamtes the noise in a non-centroided spectrum (\code{Spectrum}
      instance).
      \code{method} could be \code{"MAD"} or \code{"SuperSmoother"}.
      The arguments \code{\dots} are passed to the noise estimator functions
      implemented in \code{MALDIquant::estimateNoise}.
      Currenlty only the \code{method = "SuperSmoother"} accepts additional
      arguments, e.g. \code{span}. Please see \code{\link[stats]{supsmu}} for
      details.
      This method returns a two-column matrix with the m/z and intensity values
      in the first and the second column.
    }

    \item{\code{signature(object = "MSnExp", method = "character", ...)}}{
      Estimates noise for all spectra in \code{object}.}
  }
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  \code{\link{pickPeaks}}, and the underlying method in \code{MALDIquant}:
  \code{estimateNoise}.
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
estimateNoise(sp1, method = "SuperSmoother")
}

\keyword{methods}
