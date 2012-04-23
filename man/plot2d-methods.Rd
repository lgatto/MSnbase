\name{plot2d-methods}
\docType{methods}
\alias{plot2d-methods}
\alias{plot2d,MSnExp-method}
\alias{plot2d,data.frame-method}
\alias{plot2d}


\title{The 'plot2d' method for 'MSnExp' quality assessment}


\description{
  These methods plot the retention time vs. precursor MZ for the whole
  \code{"\linkS4class{MSnExp}"} experiment. Individual dots will be
  colour-coded to describe individual spectra's peaks count, total ion
  count, precursor charge (MS2 only) or file of origin.

  The methods make use the \code{ggplot2} system. An object of class
  'ggplot' is returned invisibly. 
}

\arguments{
  \item{object}{
    An object of class \code{"\linkS4class{MSnExp}"} or a
    \code{data.frame}. In the latter case, the data frame must have numerical
    columns named 'retention.time' and 'precursor.mz' and one of 'tic',
    'file', 'peaks.count' or 'charge', depending on the \code{z}
    parameter. Such a data frame is typically generated using the
    \code{header} method on \code{"\linkS4class{MSnExp}"} object.
  }
  \item{z}{ A character indicating according to what variable to colour
    the dots. One of, possibly abreviated, "ionCount" (total ion count),
    "file" (raw data file), "peaks.count"
    (peaks count) or "charge" (precursor charge). 
  }
  \item{alpha}{ Numeric [0,1] indicating transparence level of points. } 
  \item{plot}{ A logical indicating whether the plot should be printed
    (default is 'TRUE'). }
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", ...)}}{
      Plots a 'MSnExp' summary.
    }
    \item{\code{signature(object = "data.frame", ...)}}{
      Plots a summary of the 'MSnExp' experiment described by the data
      frame. 
    }
  }
}

\seealso{
  The \code{\link{plotDensity}} and \code{\link{plotMzDelta}} methods
  for other QC plots.
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
itraqdata
plot2d(itraqdata,z="ionCount")
plot2d(itraqdata,z="peaks.count")
plot2d(itraqdata,z="charge")
}


\keyword{methods}
