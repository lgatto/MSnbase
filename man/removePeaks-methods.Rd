\name{removePeaks-methods}
\docType{methods}
\alias{removePeaks-methods}
% \alias{removePeaks,MSnExp-method}
% \alias{removePeaks,Spectrum-method}
\alias{removePeaks}


\title{ Removes low intensity peaks }


\description{
  This method sets low intensity peaks from individual spectra
  (\code{Spectrum} instances) or whole experiments (\code{MSnExp}
  instances) to 0. The intensity threshold is set with the \code{t}
  parameter. Default is the \code{"min"} character. The threshold is
  then set as the non-0 minimum intensity found in the spectrum. Any
  other numeric values is valid. All peaks with maximum intensity
  smaller or equal to \code{t} are set to 0.

  Note that the number of peaks is not changed; the peaks below the
  threshold are set to 0 and the object is not cleanded out (see
  \code{\link{clean}}). An illustrative example is shown below. 
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", t, verbose = "logical" )}}{
      Removes low intensity peaks of all spectra in \code{MSnExp}
      object. \code{t} sets the minimum peak intensity. Default is
      "min", i.e the smallest intensity in each spectrum. Other
      \code{numeric} values are valid.
      Displays a control bar if verbose set to \code{TRUE}
      (default). Returns a new \code{MSnExp} instance. 
    }
    \item{\code{signature(object = "Spectrum", t)}}{
      Removes low intensity peaks of \code{Spectrum} object.  \code{t}
      sets the minimum peak intensity. Default is "min", i.e the
      smallest intensity in each spectrum. Other \code{numeric} values
      are valid.
      Returns a new \code{Spectrum} instance.
    }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\link{clean}} and \code{\link{trimMz}} for other spectra
  processing methods.
}

\examples{
int <- c(2,0,0,0,1,5,1,0,0,1,3,1,0,0,1,4,2,1)
sp1 <- new("Spectrum2",
           intensity=int,
           mz=1:length(int))
sp2 <- removePeaks(sp1) ## no peaks are removed here
                        ## as min intensity is 1 and
                        ## no peak has a max int <= 1
sp3 <- removePeaks(sp1,3)
intensity(sp1)
intensity(sp2)
intensity(sp3)

peaksCount(sp1) == peaksCount(sp2)
peaksCount(sp3) <= peaksCount(sp1)

data(itraqdata)
itraqdata2 <- removePeaks(itraqdata, t = 2.5e5)
table(unlist(intensity(itraqdata)) == 0)
table(unlist(intensity(itraqdata2)) == 0)
processingData(itraqdata2)
}

\keyword{methods}
