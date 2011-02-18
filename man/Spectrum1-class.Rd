\name{Spectrum1-class}
\Rdversion{1.1}
\docType{class}
\alias{Spectrum1-class}
\alias{class:Spectrum1}
\alias{Spectrum1}

\alias{polarity,Spectrum-method}
\alias{polarity}

\title{The "Spectrum1" Class for MS1 Spectra }
\description{
  \code{Spectrum1} extends the \code{"\linkS4class{Spectrum}"} class and
  introduces an MS1 specific attribute in addition to the slots in
  \code{"\linkS4class{Spectrum}"}. \code{Spectrum1} instances are not
  created directly but are contained in the \code{assayData} slot of an
  \code{"\linkS4class{MSnExp}"}. 
}

\section{Slots}{
  \describe{
    \item{\code{polarity}:}{Object of class \code{"integer"} indicating
      the polarity if the ion. }
    %% -------------------------------------------------------------- %%
    \item{\code{msLevel}:}{Object of class \code{"integer"} indicating
      the MS level: always 1 in this case (inherited from
      \code{"\linkS4class{Spectrum}"}). }
    \item{\code{peaksCount}:}{Object of class \code{"integer"}
      indicating the number of MZ peaks (inherited from
      \code{"\linkS4class{Spectrum}"}). }
    \item{\code{rt}:}{Object of class \code{"numeric"} indicating the
      retention time (in seconds) for the current ions (inherited from
      \code{"\linkS4class{Spectrum}"}). }
    \item{\code{acquisitionNum}:}{Object of class \code{"integer"}
      corresponding to the acquisition number of the current spectrum
      (inherited from \code{"\linkS4class{Spectrum}"}). }
    \item{\code{scanIndex}:}{Object of class \code{"integer"} indicating
      the scan index of the current spectrum (inherited from
      \code{"\linkS4class{Spectrum}"}). }
    \item{\code{mz}:}{Object of class \code{"numeric"} of length equal
      to the peaks count (see \code{peaksCount} slot) indicating the MZ
      values that have been measured for the current ion (inherited from
      \code{"\linkS4class{Spectrum}"}). }
    \item{\code{intensity}:}{Object of class \code{"numeric"} of same
      length as \code{mz} indicating the intensity at which each \code{mz}
      datum has been measured \code{"\linkS4class{Spectrum}"}).}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      indicating the versions of the \code{Spectrum} and \code{Spectrum1}
      classes of the instance. Intended for developer use and debugging. }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{Spectrum}"}, directly.
  Class \code{"\linkS4class{Versioned}"}, by class "Spectrum", distance 2.
}

\section{Methods}{
  See \code{"\linkS4class{Spectrum}"} for additional accessors and
  methods to process \code{Spectrum1} objects.
  \describe{
    \item{\code{polarity(object)}}{ Returns the polarity of the spectrum
      as an integer. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  Virtual super-class \code{"\linkS4class{Spectrum}"},
  \code{"\linkS4class{Spectrum2}"} for MS2 spectra and
  \code{"\linkS4class{MSnExp}"} for a full experiment container. 
}

\keyword{classes}
