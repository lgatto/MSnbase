\name{removeNoId-methods}
\docType{methods}
\alias{removeNoId-methods}
\alias{removeNoId}

\title{ Removes non-identified features }

\description{ The method removes non-identifed features in \code{MSnExp}
  and \code{MSnSet} instances using relevant information from the
  \code{feaureData} slot of a user-provide filtering vector of logicals.

}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", fcol = "pepseq", keep =
	NULL)}}{ Removes the feature from \code{object} that have a
	feature \code{fcol} (default is \code{"pepseq"}) equal to
	\code{NA}. Alternatively, one can also manually define
	\code{keep}, a vector of logical, defining the feature to be
	retained.
    }
    \item{\code{signature(object = "MSnSet", fcol = "pepseq", keep =
	NULL)}}{ As above of \code{MSnSet} instances. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\linkS4class{MSnExp}} and \code{\link{MSnSet}}.
}

\examples{
  quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "mzXML$")
  identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                     full.name = TRUE, pattern = "dummyiTRAQ.mzid")
  msexp <- readMSData(quantFile)
  msexp <- addIdentificationData(msexp, identFile)
  fData(msexp)$sequence
  length(msexp)

  ## using default fcol
  msexp2 <- removeNoId(msexp)
  length(msexp2)
  fData(msexp2)$sequence

  ## using keep
  print(fvarLabels(msexp))
  (k <- fData(msexp)$'MS.GF.EValue' > 75)
  k[is.na(k)] <- FALSE
  k
  msexp3 <- removeNoId(msexp, keep = k)
  length(msexp3)
  fData(msexp3)$sequence
}

\keyword{methods}
