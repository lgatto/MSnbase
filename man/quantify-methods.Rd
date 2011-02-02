\name{quantify-methods}
\docType{methods}
\alias{quantify-methods}
\alias{quantify}

\title{ Quantifies 'MSnExp' and 'Spectrum' objects }


\description{
  This method quantifies individual \code{"\linkS4class{Spectra}"}
  objects or full \code{"\linkS4class{MSnExp}"} experiments. Current
  implementation quantify specific peaks, defined by the
  \code{reporters} parameter using one of the possible \code{method}
  methods, which essentially qualifies for MSMS quantitation using
  isobaric tags, altough any peak can be defined in \code{reporters}.

  Eventually, each of the peaks defined in the \code{reporters} argument
  as an \code{"\linkS4class{ReporterIons}"} object will be
  quantified. This will be done for every spectrum in an experiment if
  an \code{"\linkS4class{MSnExp}"} is supplied. 

  Mention warnings.
  
}

\arguments{
  \item{object}{ Objects of class \code{"\linkS4class{Spectra}"} or
    \code{"\linkS4class{MSnExp}"}. 
  }
  \item{method}{
    Peak quantitation method. One of, possibly abreviated
    "trapezoidation", "max", or "sum". These methods return respectively
    the area under the peak, the maximum of the peak or the sum of all
    intensities of the peaks. 
  }
  \item{reporters}{
    An object of class \code{"\linkS4class{ReporterIons}"} that defines
    the peaks to be quantified.
  }
  \item{verbose}{
    Verbose of the output (only for \code{MSnExp} objects).
  }
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", method = "character",
	reporters = "ReporterIons", verbose = "logical")}}{
      Describe method and return value.
    }

    \item{\code{signature(object = "Spectrum", method = "character",
	reporters = "ReporterIons")}}{
      Describe method and return value.
    }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
## Give an example here
}


\keyword{methods}

