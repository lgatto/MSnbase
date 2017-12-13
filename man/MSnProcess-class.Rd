\name{MSnProcess-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnProcess-class}
\alias{class:MSnProcess}
\alias{MSnProcess}

\alias{fileNames,MSnProcess-method}
\alias{show,MSnProcess-method}
\alias{combine,MSnProcess,MSnProcess-method}
\alias{coerce,MSnProcess,list-method}

\title{The "MSnProcess" Class }

\description{
  \code{MSnProcess} is a container for MSnExp and MSnSet processing
  information. It records data files, processing steps, thresholds,
  analysis methods and times that have been applied to MSnExp or MSnSet
  instances. 
}

\section{Slots}{
  \describe{
    \item{\code{files}:}{Object of class \code{"character"} storing the
      raw data files used in experiment described by the
      \code{"MSnProcess"} instance. }
    \item{\code{processing}:}{Object of class \code{"character"} storing
      all the processing steps and times. }
    \item{\code{merged}:}{Object of class \code{"logical"} indicating
      whether spectra have been merged. }
    \item{\code{cleaned}:}{Object of class \code{"logical"} indicating
      whether spectra have been cleaned. See \code{\link{clean}} for
      more details and examples. } 
    \item{\code{removedPeaks}:}{Object of class \code{"character"}
      describing whether peaks have been removed and which threshold was
      used. See \code{\link{removePeaks}} for more details and examples. }
    \item{\code{smoothed}:}{Object of class \code{"logical"} indicating
      whether spectra have been smoothed. }
    \item{\code{trimmed}:}{Object of class \code{"numeric"} documenting
      if/how the data has been trimmed. }
    \item{\code{normalised}:}{Object of class \code{"logical"}
      describing whether and how data have been normalised. }
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
  \describe{
    \item{fileNames}{\code{signature(object = "MSnProcess")}: Returns the
      file names used in experiment described by the \code{"MSnProcess"}
      instance. } 
    \item{show}{\code{signature(object = "MSnProcess")}: Displays object
      content as text. }
    \item{combine}{\code{signature(x = "MSnProcess", y = "MSnProcess")}:
      Combines multiple \code{MSnProcess} instances. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\note{
  This class is likely to be updated using an \code{AnnotatedDataFrame}.
}


\seealso{
  See the \code{"\linkS4class{MSnExp}"} and \code{"\linkS4class{MSnSet}"}
  classes that actually use \code{MSnProcess} as a slot. 
}

\examples{
showClass("MSnProcess")
}

\keyword{classes}
