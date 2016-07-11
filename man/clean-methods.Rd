\name{clean-methods}
\docType{methods}
\alias{clean-methods}
% \alias{clean,MSnExp-method}
% \alias{clean,Spectrum-method}
\alias{clean}


\title{ Cleans 'MSnExp' or 'Spectrum' instances }

\description{
  This method cleans out individual spectra (\code{Spectrum} instances)
  or whole experiments (\code{MSnExp} instances) of 0-intensity
  peaks. Unless \code{all} is set to \code{FALSE}, origninal 0-intensity
  values are retained only around peaks. If more than two 0's were
  separating two peaks, only the first and last ones, those directly
  adjacent to the peak ranges are kept. If two peaks are separated by
  only one 0-intensity value, it is retained. An illustrative example is
  shown below.  
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", all = "logical" verbose =
	"logical")}}{ Cleans all spectra in \code{MSnExp}
      object. Displays a control bar if verbose set to \code{TRUE}
      (default). Returns a cleaned \code{MSnExp} instance.
    }

    \item{\code{signature(object = "Spectrum", all = "logical",
	msLevel. = "numeric")}}{ Cleans the \code{Spectrum}
	object. Returns a cleaned \code{Spectrum} instance. If \code{all
	= TRUE}, then all zeros are removed. \code{msLevel.} defines the
	level of the spectrum, and if \code{msLevel(object) !=
	msLevel.}, cleaning is ignored. Only relevant when called from
	\code{OnDiskMSnExp} and is only relevant for developers.  }
    
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\link{removePeaks}} and \code{\link{trimMz}} for other spectra
  processing methods.
}

\examples{
int <- c(1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0)
sp1 <- new("Spectrum2",
           intensity=int,
           mz=1:length(int))
sp2 <- clean(sp1) ## default is all=FALSE
intensity(sp1)
intensity(sp2)
intensity(clean(sp1, all = TRUE))

mz(sp1)
mz(sp2)
mz(clean(sp1, all = TRUE))

data(itraqdata)
itraqdata2 <- clean(itraqdata)
sum(peaksCount(itraqdata))
sum(peaksCount(itraqdata2))
processingData(itraqdata2)
}


\keyword{methods}
