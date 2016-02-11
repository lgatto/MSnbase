\name{Spectrum1-class}
\Rdversion{1.1}
\docType{class}
\alias{Spectrum1-class}
\alias{class:Spectrum1}
\alias{Spectrum1}

\title{The "Spectrum1" Class for MS1 Spectra }
\description{
  \code{Spectrum1} extends the \code{"\linkS4class{Spectrum}"} class and
  introduces an MS1 specific attribute in addition to the slots in
  \code{"\linkS4class{Spectrum}"}. \code{Spectrum1} instances are not
  created directly but are contained in the \code{assayData} slot of an
  \code{"\linkS4class{MSnExp}"}. 
}

\section{Slots}{
  See the \code{"\linkS4class{Spectrum}"} class for inherited slots.
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
