\name{extractSpectra-methods}
\docType{methods}
\alias{extractSpectra-methods}
\alias{extractSpectra}


\title{ Extracts a set of spectra from an 'MSnExp' object }

\description{
  Extracts the spectra defined by the \code{selected} parameter and
  returns a subsetted \code{"\linkS4class{MSnExp}"} instance. The
  \code{featureData} slot is also subsetted accordingly.

  Note that values in \code{selected} are recycled shorter than
  \code{length(MSnExp)}. 
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", selected = "logical")}}{
      Returns an \code{"\linkS4class{MSnExp}"} object containing spectra 
      defined by the logical \code{selected} argument.
    }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\link{extractPrecSpectra}}
}

\examples{
file <- dir(system.file(package="MSnbase",dir="extdata"),
            full.name=TRUE,pattern="mzXML$")
aa <- readMzXMLData(file,verbose=FALSE)
sel <- rep(FALSE,length(aa))
sel[1:3] <- TRUE
bb <- extractSpectra(aa,sel)
stopifnot(length(bb) == sum(sel))
processingData(bb)
}

\keyword{methods}
