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
  For noisy spectra there are currently two different noise estimators
  available, the Median Absolute Deviation (\code{method = "MAD"}) and
  Friedman's Super Smoother (\code{method = "SuperSmoother"}),
  as implemented in the \code{MALDIquant::detectPeaks} and
  \code{MALDIquant::estimateNoise} functions respectively.

  The method supports also to optionally \emph{refine} the m/z value of
  the identified centroids by considering data points that belong (most
  likely) to the same mass peak. The m/z value is calculated as an
  intensity weighted average of the m/z values within the peak region. 
  How the peak region is defined depends on the method chosen:

  \code{refineMz = "kNeighbors"}: m/z values (and their respective
  intensities) of the\code{2 * k} closest signals to the centroid are
  used in the intensity weighted average calculation. The number of
  neighboring signals can be defined with the argument \code{k}.

  \code{refineMz = "descendPeak"}: the peak region is defined by
  descending from the identified centroid/peak on both sides until the
  measured signal increases again. Within this defined region all
  measurements with an intensity of at least \code{signalPercentage} of
  the centroid's intensity are used to calculate the refined m/z. By
  default the descend is stopped when the first signal that is equal or
  larger than the last observed one is encountered. Setting
  \code{stopAtTwo = TRUE}, two consecutively increasing signals are
  required.
  
  By default (\code{refineMz = "none"}, simply the m/z of the largest
  signal (the identified centroid) is reported. See below for examples.
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "MSnExp", halfWindowSize = "integer",
        method = "character", SNR = "numeric", verbose = "logical",
	refineMz = "character", ...)}}{
      Performs the peak picking for all spectra in an \code{MSnExp} instance.
      \code{method} could be \code{"MAD"} or \code{"SuperSmoother"}.
      \code{halfWindowSize} controls the window size of the peak picking
      algorithm. The resulting window size is \code{2 * halfWindowSize + 1}.
      The size should be nearly (or slightly larger) the \emph{FWHM}
      (full width at half maximum).
      A local maximum is considered as peak if its intensity is \code{SNR}
      times larger than the estimated noise.
      \code{refineMz} allows to choose a method for an optional centroid
      m/z refinement (see description for more details). Choises are
      \code{"none"} (default, no m/z refinement), \code{"kNeighbors"}
      and \code{"descendPeak"}.
      The arguments \code{\dots} are passed to the noise estimator or
      m/z refinement functions.
      For the noise estimator functions, currenlty only the
      \code{method = "SuperSmoother"} accepts additional arguments,
      e.g.  \code{span}. Please see \code{\link[stats]{supsmu}} for
      details. \code{refineMethod = "kNeighbors"} supports additional
      argument \code{k} and \code{refineMethod = "descendPeak"}
      arguments \code{signalPercentage} and \code{stopAtTwo}. See
      description above for more details.
      
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
  Sebastian Gibb <mail@sebastiangibb.de> with contributions from
  Johannes Rainer.
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


## Examples for refineMz:
ints <- c(5, 3, 2, 3, 1, 2, 4, 6, 8, 11, 4, 7, 5, 2, 1, 0, 1, 0, 1, 1, 1, 0)
mzs <- 1:length(ints)
sp1 <- new("Spectrum1", intensity = ints, mz = mzs, centroided = FALSE)
plot(mz(sp1), intensity(sp1), type = "h")

## Without m/z refinement:
sp2 <- pickPeaks(sp1)           
points(mz(sp2), intensity(sp2), col = "darkgrey")
## Using k = 1, closest signals
sp3 <- pickPeaks(sp1, refineMz = "kNeighbors", k = 1)           
points(mz(sp3), intensity(sp3), col = "green", type = "h")

## Using descendPeak requiring at least 50% or the centroid's intensity
sp4 <- pickPeaks(sp1, refineMz = "descendPeak", signalPercentage = 50)           
points(mz(sp4), intensity(sp4), col = "red", type = "h")
           
}

\keyword{methods}
