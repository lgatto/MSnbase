\name{Spectrum-class}
\Rdversion{1.1}
\docType{class}
\alias{Spectrum-class}
\alias{class:Spectrum}
\alias{Spectrum}

%% -- accessors -- %%
\alias{acquisitionNum,Spectrum-method}
\alias{acquisitionNum}
\alias{scanIndex,Spectrum-method}
\alias{scanIndex}
\alias{fromFile,Spectrum-method}
\alias{fromFile}
\alias{intensity,Spectrum-method}
\alias{intensity}
\alias{msLevel,Spectrum-method}
\alias{msLevel}
\alias{mz,Spectrum-method}
\alias{mz}
\alias{peaksCount,Spectrum,missing-method}
\alias{peaksCount}
\alias{rtime,Spectrum-method}
\alias{rtime}
\alias{tic,Spectrum-method}
\alias{tic}
\alias{ionCount,Spectrum-method}
\alias{ionCount}
\alias{centroided,Spectrum-method}
\alias{centroided}
\alias{centroided<-,Spectrum,ANY-method}
\alias{centroided<-,Spectrum,logical-method}
\alias{centroided<-}

%% -- other -- %%
\alias{bin,Spectrum-method}
\alias{clean,Spectrum-method}
\alias{compareSpectra,Spectrum,Spectrum-method}
\alias{pickPeaks,Spectrum-method}
\alias{quantify,Spectrum-method}
\alias{quantify,Spectrum,character-method}
\alias{removePeaks,Spectrum-method}
\alias{show,Spectrum-method}
\alias{smooth,Spectrum-method}
\alias{trimMz,Spectrum-method}
\alias{trimMz,Spectrum,numeric-method}
\alias{isEmpty,Spectrum-method}

\alias{coerce,Spectrum,data.frame-method}
\alias{as.data.frame.Spectrum}

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
    \item{\code{tic}:}{Object of class \code{"numeric"} indicating the
      total ion current.}
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
    \item{\code{centroided}:}{Object of class \code{"logical"}
      indicating if instance is centroided ('TRUE') of uncentroided
      ('FALSE'). }
    \item{\code{fromFile}:}{Object of class \code{"integer"} referencing
      the file the spectrum originates. The file names are stored in the
      \code{processingData} slot of the \code{"\linkS4class{MSnExp}"} or
      \code{"\linkS4class{MSnSet}"} instance that contains the current
      \code{"Spectrum"} instance. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      indicating the version of the \code{Spectrum} class. Intended for
      developer use and debugging.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

\section{Methods}{
  \describe{
    \item{\code{acquisitionNum(object)}}{ Returns the acquisition number
      of the spectrum as an integer. }
    \item{\code{scanIndex(object)}}{ Returns the scan index of the
      spectrum as an integer. }
    \item{\code{centroided(object)}}{ Indicates whether spectrum is
      centroided ('TRUE') or uncentroided ('FALSE'). }
    \item{\code{centroided(object) <- value}}{ Sets the 'centroided'
      status of the spectrum object. }
    \item{\code{fromFile(object)}}{ Returns the index of the raw data
      file from which the current instances originates as an
      integer. }
    \item{\code{intensity(object)}}{ Returns an object of class
      \code{"numeric"} containing the intensities of the spectrum. }
    \item{\code{msLevel(object)}}{ Returns an MS level of the spectrum
      as an integer. }
    \item{\code{mz(object, ...)}}{ Returns an object of class
      \code{"numeric"} containing the MZ value of the spectrum
      peaks. Additional arguments are currently ignored. } 
    \item{\code{peaksCount(object)}}{ Returns the number of peaks
      (possibly of 0 intensity) as an integer. }
    \item{\code{rtime(object, ...)}}{ Returns the retention time for the
      spectrum as an integer. Additional arguments are currently
      ignored. }
    \item{\code{ionCount(object)}}{ Returns the total ion count for the
      spectrum as a numeric. }
    \item{\code{tic(object, ...)}}{ Returns the total ion current for the
      spectrum as a numeric. Additional arguments are currently ignored. }
    %% -------------------------------------------------- %%
    \item{bin}{\code{signature(object = "Spectrum")}: Bins Spectrum.
      See \code{\link{bin}} documentation for more details and examples. }
    \item{clean}{\code{signature(object = "Spectrum")}: Removes unused 0
      intensity data points. See \code{\link{clean}} documentation
      for more details and examples. }
    \item{compareSpectra}{\code{signature(object1 = "Spectrum",
      object2 = "Spectrum")}: Compares spectra.  See
      \code{\link{compareSpectra}} documentation for more details and
      examples. }
    \item{pickPeaks}{\code{signature(object = "Spectrum")}: Performs the peak
      picking to generate a centroided spectrum.
      See \code{\link{pickPeaks}} documentation for more
      details and examples. }
    \item{plot}{\code{signature(x = "Spectrum", y = "missing")}: Plots
      intensity against mz.
      See \code{\link{plot.Spectrum}} documentation for more details. }
    \item{plot}{\code{signature(x = "Spectrum", y = "Spectrum")}: Plots
      two spectra above/below each other.
      See \code{\link{plot.Spectrum.Spectrum}} documentation for more
      details.}
    \item{plot}{\code{signature(x = "Spectrum", y = "character")}: Plots
      an MS2 level spectrum and its highlight the fragmention peaks.
      See \code{\link{plot.Spectrum.character}} documentation for more
      details.}
    \item{quantify}{\code{signature(object = "Spectrum")}: Quatifies
      defined peaks in the spectrum.
      See \code{\link{quantify}} documentation for more details. }
    \item{removePeaks}{\code{signature(object = "Spectrum")}: Remove
      peaks lower that a threshold \code{t}. See
      \code{\link{removePeaks}} documentation for more details and
      examples. }
    \item{smooth}{\code{signature(x = "Spectrum")}: Smooths spectrum.
      See \code{\link{smooth}} documentation for more details and examples. }
    \item{show}{\code{signature(object = "Spectrum")}: Displays object
      content as text. }
    \item{trimMz}{\code{signature(object = "Spectrum")}: Trims the MZ
      range of all the spectra of the \code{MSnExp} instance. See
      \code{\link{trimMz}} documentation for more details and
      examples. }
    \item{isEmpty}{\code{signature(x = "Spectrum")}: Checks if the
      \code{x} is an empty \code{Spectrum}. }   
    %% -------------------------------------------------- %%
    \item{as}{\code{signature(object = "Spectrum", "data.frame")}:
      Coerces the \code{Spectrum} object to a two-column
      \code{data.frame} containing intensities and MZ values. }
  }
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
