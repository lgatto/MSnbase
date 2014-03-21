\name{addIdentificationData-methods}
\docType{methods}
\alias{addIdentificationData-methods}
\alias{addIdentificationData}

\title{ Adds Identification Data }

\description{
  This methods add identification data to an experiment
  \code{"\linkS4class{MSnExp}"} or to a \code{"\linkS4class{MSnSet}"}.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", filenames = "character",
	    verbose = "logical")}}{ 
      Adds the identification data stored in mzIdentML files to a 
      \code{"\linkS4class{MSnExp}"} instance. The method handles
      one or multiple mzIdentML files provided via \code{filenames}.
      The \code{verbose} argument (default is
      \code{TRUE}) defines whether status messages should be showed.
    }
    \item{\code{signature(object = "MSnSet", filenames = "character",
	    verbose = "logical")}}{ 
      Adds the identification data stored in mzIdentML files to an
      \code{"\linkS4class{MSnSet}"} instance. Please note that that would be
      only possible if the \code{"\linkS4class{MSnSet}"} was generated from an
      \code{"\linkS4class{MSnExp}"} object or the \code{data.frame} in 
      the \code{featureData} slot has the columns \code{file} and
      \code{acquisition.number}.\cr
      The method handles one or multiple mzIdentML files provided via 
      \code{filenames}. The \code{verbose} argument (default is
      \code{TRUE}) defines whether status messages should be showed.
    }
  }
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  \code{\linkS4class{MSnExp}} and \code{\link{MSnSet}}.
}

\examples{
## find path to a mzXML file
quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
## find path to a mzIdentML file
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzid$")

## create basic MSnExp
msexp <- readMSData(quantFile)

## add identification information
msexp <- addIdentificationData(msexp, identFile)
}

\keyword{methods}
