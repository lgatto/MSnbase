\name{removeReporters-methods}
\docType{methods}
\alias{removeReporters-methods}
\alias{removeReporters}

\title{ Removes reporter ion tag peaks }


\description{
  This methods sets all the reporter tag ion peaks from one MS2
  spectrum or all the MS2 spectra of an experiment to 0. Reporter data
  is specified using an \code{"\linkS4class{ReporterIons}"}
  instance. The peaks are selected around the expected reporter ion
  m/z value +/- the reporter width. 
  Optionally, the spectrum/spectra can be \code{cleaned} to
  remove successive 0 intensity data points (see the \code{\link{clean}}
  function for details).

  Note that this method only works for MS2 spectra or experiments that
  contain MS2 spectra. It will fail for MS1 spectrum.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", reporters = "ReporterIons",
	clean = "logical", verbose = "logical" )}}{ The reporter ion
      peaks defined in the \code{reporters} instance of all the MS2
      spectra of the \code{"\linkS4class{MSnExp}"} instance are set to 0
      and, if \code{clean} is set to \code{TRUE}, cleaned. The default
      value of \code{reporters} is \code{NULL}, which leaves the spectra
      as unchanged. The \code{verbose} parameter (default is
      \code{TRUE}) defines whether a progress bar should be showed.
    }
    \item{\code{signature(object = "Spectrum", reporters = "ReporterIons",
	clean = "FALSE")}}{ The reporter ion peaks defined in the
      \code{reporters} instance of MS2 \code{"\linkS4class{Spectrum}"}
      instance are set to 0 and, if \code{clean} is set to \code{TRUE},
      cleaned. The default value of \code{reporters} is \code{NULL},
      which leaves the spectrum as unchanged.
    }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\link{clean}} and \code{\link{removePeaks}} for other spectra
  processing methods.
}

\examples{
sp1 <- itraqdata[[1]]
sp2 <- removeReporters(sp1,reporters=iTRAQ4)
sel <- mz(sp1) > 114 & mz(sp1) < 114.2
mz(sp1)[sel]
intensity(sp1)[sel]
plot(sp1,full=TRUE,reporters=iTRAQ4)
intensity(sp2)[sel]
plot(sp2,full=TRUE,reporters=iTRAQ4)
}

\keyword{methods}
