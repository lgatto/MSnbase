\name{Spectrum2-class}
\Rdversion{1.1}
\docType{class}
\alias{Spectrum2-class}
\alias{class:Spectrum2}
\alias{Spectrum2}

\alias{precursorMz,Spectrum-method}
\alias{precursorMz}
\alias{precScanNum,Spectrum-method}
\alias{precScanNum}
\alias{precAcquisitionNum,Spectrum-method}
\alias{precAcquisitionNum}
\alias{precursorIntensity,Spectrum-method}
\alias{precursorIntensity}
\alias{precursorCharge,Spectrum-method}
\alias{precursorCharge}
\alias{collisionEnergy,Spectrum-method}
\alias{collisionEnergy}
\alias{removeReporters,Spectrum-method}
\alias{calculateFragments,character,Spectrum2-method}

\title{The "Spectrum2" Class for MSn Spectra}

\description{
  
  \code{Spectrum2} extends the \code{"\linkS4class{Spectrum}"} class and
  introduces several MS2 specific attributes in addition to the slots in
  \code{"\linkS4class{Spectrum}"}. Since version 1.99.2, this class is
  used for any MS levels > 1. \code{Spectrum2} are not created directly
  but are contained in the \code{assayData} slot of an
  \code{"\linkS4class{MSnExp}"}.

  In version 1.19.12, the \code{polarity} slot had been added to the
  \code{"\linkS4class{Spectrum}"} class (previously in
  \code{"\linkS4class{Spectrum1}"}). Hence, \code{"Spectrum2"} objects
  created prior to this change will not be valid anymore, since they
  will miss the \code{polarity} slots. Object can be appropriately
  updated using the \code{updateObject} method.
  
}


\section{Slots}{
  See the \code{"\linkS4class{Spectrum}"} class for inherited slots.
  \describe{
    \item{\code{merged}:}{Object of class \code{"numeric"} indicating of
      how many combination the current spectrum is the result of. }
    \item{\code{precScanNum}:}{Object of class \code{"integer"} indicating
      the precursor MS scan index in the original input file. Accessed
      with the \code{precScanNum} or \code{precAcquisitionNum} methods.}
    \item{\code{precursorMz}:}{Object of class \code{"numeric"}
      providing the precursor ion MZ value. }
    \item{\code{precursorIntensity}:}{Object of class \code{"numeric"}
      providing the precursor ion intensity. }
    \item{\code{precursorCharge}:}{Object of class \code{"integer"}
      indicating the precursor ion charge. }
    \item{\code{collisionEnergy}:}{Object of class \code{"numeric"}
      indicating the collision energy used to fragment the parent ion. }
  }
}

\section{Extends}{
  Class \code{"\linkS4class{Spectrum}"}, directly.
  Class \code{"\linkS4class{Versioned}"}, by class "Spectrum", distance 2.
}

\section{Methods}{
  See \code{"\linkS4class{Spectrum}"} for additional accessors and
  methods for \code{Spectrum2} objects.
  \describe{
    \item{\code{precursorMz(object)}}{ Returns the precursor MZ value as
      a numeric. }
    \item{\code{precursorMz(object)}}{ Returns the precursor scan number
      in the original data file as an integer. }
    \item{\code{precursorIntensity(object)}}{ Returns the precursor
      intensity as a numeric. }
    \item{\code{precursorCharge(object)}}{ Returns the precursor
      intensity as a integer. }
    \item{\code{collisionEnergy(object)}}{ Returns the collision energy
      as an numeric. }
    \item{\code{removeReporters(object, ...)}}{ Removes all reporter ion
      peaks. See \code{\link{removeReporters}} documentation for more
      details and examples. }
    \item{\code{precAcquisitionNum}:}{Returns the precursor's
      acquisition number. }
    \item{\code{precScanNum}:}{See \code{precAcquisitionNum}. }
    \item{calculateFragments}{\code{signature(sequence = "character",
      object = "Spectrum2")}:
      Calculates and matches the theoretical fragments of a peptide
      \code{sequence} with the ones observed in a spectrum.
      See \code{\link{calculateFragments}} documentation
      for more details and examples. }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  Virtual super-class \code{"\linkS4class{Spectrum}"},
  \code{"\linkS4class{Spectrum1}"} for MS1 spectra and
  \code{"\linkS4class{MSnExp}"} for a full experiment container.
}

\keyword{classes}
