\name{plotDensity-methods}
\docType{methods}
\alias{plotDensity-methods}
\alias{plotDensity,MSnExp-method}
\alias{plotDensity,data.frame-method}
\alias{plotDensity}


\title{The 'plotDensity' method for 'MSnExp' quality assessment}

\description{
  These methods plot the distribution of several parameters of interest
  for the different precursor charges for \code{"\linkS4class{MSnExp}"}
  experiment.
  
  The methods make use the \code{ggplot2} system. An object of class
  'ggplot' is returned invisibly. 
}

\arguments{
  \item{object}{
    An object of class \code{"\linkS4class{MSnExp}"} or and
    'data.frame'. In the latter case, the data frame must have numerical
    columns named 'charge' and one of 'precursor.mz', 'peaks.count' or 'ionCount',
    depending on the \code{z} parameter.
    Such a data frame is typically generated using the
    \code{header} method on \code{"\linkS4class{MSnExp}"} object.
  }
  \item{z}{ A character indicating which parameter's densitiy to plot.
    One of, possibly abreviated, "ionCount" (total ion count),
    "peaks.count" (peaks count) or "precursor.mz" (precursor MZ). 
  }
  \item{log}{ Logical, whether to log transform the data (default is
    'FALSE'). } 
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
  The \code{\link{plot2d}} and \code{\link{plotDensity}} methods for
  other QC plots.
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
itraqdata
plotDensity(itraqdata,z="ionCount")
plotDensity(itraqdata,z="peaks.count")
plotDensity(itraqdata,z="precursor.mz")
}


\keyword{methods}
