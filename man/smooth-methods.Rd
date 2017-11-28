\name{smooth-methods}
\docType{methods}
\alias{smooth-methods}
% \alias{smooth,MSnExp-method}
% \alias{smooth,Spectrum-method}
\alias{smooth}


\title{ Smooths 'MSnExp' or 'Spectrum' instances }

\description{
  This method smooths individual spectra (\code{Spectrum} instances)
  or whole experiments (\code{MSnExp} instances).
  Currently, the Savitzky-Golay-Smoothing (\code{method = "SavitzkyGolay"})
  and the Moving-Average-Smoothing (\code{method = "MovingAverage"}) are 
  available, as implemented in the \code{MALDIquant::smoothIntensity} function. 
  Additional methods might be added at a later stage.
}

\section{Methods}{

  \describe{

    \item{\code{signature(x = "MSnExp", method = "character",
	halfWindowSize = "integer", verbose = "logical", ...)}}{

      Smooths all spectra in \code{MSnExp}. \code{method} could be
      \code{"SavitzkyGolay"} or
      \code{"MovingAverage"}. \code{"halfWindowSize"} controls the
      window size of the filter. The resulting window size is \code{2 *
      halfWindowSize + 1}. The best size differs depending on the
      selected \code{method}. For \code{method = "SavitzkyGolay"} it
      should be lower than \emph{FWHM} of the peaks (full width at half
      maximum; please find details in Bromba and Ziegler 1981).  The
      arguments \code{\dots} are passed to the internal functions.
      For \code{method="MovingAverage"} there is an additional \code{weighted}
      argument (default: \code{FALSE}) to indicate if the average should
      be equal weight (default) or if it should have weights depending
      on the distance from the center as calculated as
      \code{1/2^abs(-halfWindowSize:halfWindowSize)} with the sum
      of all weigths normalized to 1.
      For \code{method="SavitzkyGolay"} an additonal argument
      is \code{polynomialOrder} (default: 3). It controls the
      polynomial order of the Savitzky-Golay Filter.
      This method displays a progress bar if \code{verbose = TRUE}.
      Returns an \code{MSnExp} instance with smoothed spectra.
      }

      \item{\code{signature(x = "Spectrum", method = "character",
	  halfWindowSize = "integer", ...)}}{

	Smooths the spectrum (\code{Spectrum} instance). This method is
	the same as above but returns a smoothed \code{Spectrum} instead
	of an \code{MSnExp} object. It has no \code{verbose}
	argument. Please read the details for the above \code{MSnExp}
	method.

      } 

    }

  }

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  \code{\link{clean}}, \code{\link{pickPeaks}}, \code{\link{removePeaks}} and 
  \code{\link{trimMz}} for other spectra processing methods.
}

\references{
A. Savitzky and M. J. Golay. 1964.
Smoothing and differentiation of data by simplified least squares procedures.
Analytical chemistry, 36(8), 1627-1639.

M. U. Bromba and H. Ziegler. 1981.
Application hints for Savitzky-Golay digital smoothing filters.
Analytical Chemistry, 53(11), 1583-1586.

S. Gibb and K. Strimmer. 2012.
MALDIquant: a versatile R package for the analysis of mass spectrometry data.
Bioinformatics 28: 2270-2271.
\url{http://strimmerlab.org/software/maldiquant/}
}

\examples{
sp1 <- new("Spectrum1",
           intensity = c(1:6, 5:1),
           mz = 1:11)
sp2 <- smooth(sp1, method = "MovingAverage", halfWindowSize = 2)
intensity(sp2)

data(itraqdata)
itraqdata2 <- smooth(itraqdata, 
                     method = "MovingAverage", 
                     halfWindowSize = 2)
processingData(itraqdata2)
}

\keyword{methods}
