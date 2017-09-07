\name{trimMz-methods}
\docType{methods}
\alias{trimMz-methods}
% \alias{trimMz,MSnExp-method}
% \alias{trimMz,Spectrum-method}
\alias{trimMz}
\alias{filterMz}
\alias{filterMz-methods}


\title{ Trims 'MSnExp' or 'Spectrum' instances }


\description{
  
  This method selects a range of MZ values in a single spectrum
  (\code{Spectrum} instances) or all the spectra of an experiment
  (\code{MSnExp} instances). The regions to trim are defined by the
  range of \code{mz} argument, such that MZ values <= \code{min(mz)} and
  MZ values >= \code{max(mz)} are trimmed away.
  
}
\section{Methods}{
  \describe{
    
    \item{\code{signature(object = "MSnExp", mz = "numeric", msLevel. =
	"numeric")}}{ Trims all spectra in \code{MSnExp} object according
      to \code{mz}. If \code{msLevel.} is defined, then only spectra of
      that level are trimmer. }

    \item{\code{signature(object = "Spectrum", mz = "numeric",
	msLevel. = "numeric")}}{Trims the \code{Spectrum} object and
      retruns a new trimmed object. \code{msLevel.} defines the level
      of the spectrum, and if \code{msLevel(object) != msLevel.},
      cleaning is ignored. Only relevant when called from
      \code{OnDiskMSnExp} and is only relevant for developers.  }

  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\link{removePeaks}} and \code{\link{clean}} for other spectra
  processing methods.
}

\examples{
mz <- 1:100
sp1 <- new("Spectrum2",
           mz = mz,
           intensity = abs(rnorm(length(mz))))

sp2 <- trimMz(sp1, c(25, 75))
range(mz(sp1))
range(mz(sp2))

data(itraqdata)
itraqdata2 <- filterMz(itraqdata, c(113, 117))
range(mz(itraqdata))
range(mz(itraqdata2))
processingData(itraqdata2)
}


\keyword{methods}

