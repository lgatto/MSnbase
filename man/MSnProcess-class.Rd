\name{MSnProcess-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnProcess-class}
\alias{MSnProcess}

\title{The "MSnProcess" Class }

\description{
  \code{MSnProcess} is a container for MSnExp and MSnSet processing
  information. It records processing steps, thresholds, analysis methods
  and times that hav been applied to MSnExp or MSnSet instances.
}

\section{Slots}{
  \describe{
    \item{\code{processing}:}{Object of class \code{"character"} storing
      all the processing steps and times. }
    \item{\code{merged}:}{Object of class \code{"logical"} indicating
      whether spectra have been merged. }
    \item{\code{cleaned}:}{Object of class \code{"logical"} indicating
      whether spectra have been cleaned. }
    \item{\code{removedPeaks}:}{Object of class \code{"character"}
      describing whether peaks have been removed and which threshold was
      used. }
    \item{\code{smoothed}:}{Object of class \code{"logical"} indicating
      whether spectra have been smoothed. }
    \item{\code{centroided}:}{Object of class \code{"logical"}
      indicating whether spectra have been centroided. }
    \item{\code{normalised}:}{Object of class \code{"logical"}
      describing whether and how data have been normalised. }
    \item{\code{xcmsVersion}:}{Object of class \code{"character"}
      indicating the version of xcms used when data was loaded. }
    \item{\code{MSnbaseVersion}:}{Object of class \code{"character"}
      indicating the version of MSnbase. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      indicating the version of the \code{MSnProcess} instance. Intended for
      developer use and debugging. }
      
  }
}
\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}
\section{Methods}{
  No methods defined with class "MSnProcess" in the signature.
}
\author{
  Laurent Gatto <lg390@cam.ac.uk>
}
\note{
  This class is likely to be updated.
}


\seealso{
  See \code{"\linkS4class{MSnExp}"} and \code{"\linkS4class{MSnSet}"}
  classes that actually use \code{MSnProcess} as a slot. 
}

\examples{
showClass("MSnProcess")
}

\keyword{classes}
