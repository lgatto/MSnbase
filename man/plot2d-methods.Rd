\name{plot2d-methods}
\docType{methods}
\alias{plot2d-methods}
\alias{plot2d,data.frame,character-method}
\alias{plot2d}


\title{ Graphical overview of the 'MSnExp' object }


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
    An object of class \code{"\linkS4class{MSnExp}"} or and
    'data.frame'. In the latter case, the data frame must have numerical
    columns named 'retention.time' and 'precursor.mz' and one of 'tic',
    'file', 'peaks.count' or 'charge', depending on the \code{z}
    parameter. Such a data.frame is typically generated using the
    \code{header} method on \code{"\linkS4class{MSnExp}"} object.
  }
  \item{z}{ A character indicating according to what variable to colour
    the dots. One of, possibly abreviated, "tic" (total ion count),
    "file" (raw data file), "peaks.count"
    (peaks count) or "charge" (precursor charge). 
  }
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", z = "character")}}{
      Plots a 'MSnExp' summary.
    }
    \item{\code{signature(object = "data.frame", z = "character")}}{
      Plots a summary of the 'MSnExp' experiment described by the data
      frame. 
    }
  }
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
aa <- readMzXMLData(file,verbose=FALSE)
plot2d(aa,z="tic")
plot2d(aa,z="file")
}


\keyword{methods}
