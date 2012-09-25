\name{trimMz-methods}
\docType{methods}
\alias{trimMz-methods}
% \alias{trimMz,MSnExp-method}
% \alias{trimMz,Spectrum-method}
\alias{trimMz}


\title{ Trims 'MSnExp' or 'Spectrum' instances }


\description{
  This method selects a range of MZ values in a single spectrum
  (\code{Spectrum} instances) or all the spectra of an experiment
  (\code{MSnExp} instances). The regions to trim are defined by the
  range of \code{mzlim} argument, such that MZ values < \code{min(mzlim)}
  and MZ values > \code{max(mzlim)} are trimmed away.
}
\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", mzlim = "numeric")}}{
      Trims all spectra in \code{MSnExp} object according to
      \code{mzlim}. Returns a cleaned \code{MSnExp} instance.
    }
    \item{\code{signature(object = "Spectrum", mzlim = "numeric")}}{
      Trims the \code{Spectrum} object and retruns a new trimmed object.
    }
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
           mz=mz,
           intensity=abs(rnorm(length(mz))))
sp2 <- trimMz(sp1,c(25,75))
range(mz(sp1))
range(mz(sp2))

data(itraqdata)
itraqdata2 <- trimMz(itraqdata,c(113,117))
range(mz(itraqdata))
range(mz(itraqdata2))
processingData(itraqdata2)
}


\keyword{methods}

