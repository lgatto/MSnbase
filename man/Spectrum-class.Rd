\name{Spectrum-class}
\Rdversion{1.1}
\docType{class}
\alias{Spectrum-class}
\alias{class:Spectrum}
\alias{Spectrum}

\title{The "Spectrum" Class}

\description{
  Virtual container for spectrum data common to all different types of
  spectra. A \code{Spectrum} object can not be directly instanciated. Use
  \code{"\linkS4class{Spectrum1}"} and \code{"\linkS4class{Spectrum2}"}
  instead.  
}

\section{Slots}{
  \describe{
    \item{\code{msLevel}:}{Object of class \code{"integer"} indicating
      the MS level: 1 for MS1 level \code{Spectrum1} objects and 2 for MSMSM
      \code{Spectrum2} objects. Levels > 2 have not been tested and will be
      handled as MS2 spectra. }
    \item{\code{peaksCount}:}{Object of class \code{"integer"}
      indicating the number of MZ peaks. }
    \item{\code{rt}:}{Object of class \code{"numeric"} indicating the
      retention time (in seconds) for the current ions. }
    \item{\code{acquisitionNum}:}{Object of class \code{"integer"}
      corresponding to the acquisition number of the current spectrum. }
    \item{\code{scanIndex}:}{Object of class \code{"integer"} indicating
      the scan index of the current spectrum. }
    \item{\code{mz}:}{Object of class \code{"numeric"} of length equal
      to the peaks count (see \code{peaksCount} slot) indicating the MZ
      values that have been measured for the current ion.  }
    \item{\code{intensity}:}{Object of class \code{"numeric"} of same
      length as \code{mz} indicating the intensity at which each \code{mz}
      datum has been measured. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      indicating the version of the \code{Spectrum} class. Intended for
      developer use and debugging.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}
\note{
  This is a virtual class and can not be instanciated directly.
}

\seealso{
  Instaciable sub-classes \code{"\linkS4class{Spectrum1}"} and
  \code{"\linkS4class{Spectrum2}"} for MS1 and MS2 spectra.
}

\keyword{classes}
