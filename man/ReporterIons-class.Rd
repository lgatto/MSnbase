\name{ReporterIons-class}
\Rdversion{1.1}
\docType{class}
\alias{ReporterIons-class}
\alias{class:ReporterIons}
\alias{ReporterIons}

\title{The "ReporterIons" Class}

\description{
  The \code{ReporterIons} class allows to define a set of isobaric
  reporter ions that are used for quantification in MSMS
  mode, e.g. iTRAQ (isobaric tag for relative and absolute quantitation)
  or TMT (tandem mass tags).
  \code{ReporterIons} instances can them be used when quantifying
  \code{"\linkS4class{MSnExp}"} data of plotting the reporters peaks
  based on in \code{"\linkS4class{Spectrum2}"} ojects.

  Some reporter ions are provided with \code{MSnbase} an can be loaded
  with the \code{\link{data}} function. These reporter ions data sets
  are:
  \describe{
    \item{\code{iTRAQ4}:}{\code{ReporterIon} object for the iTRAQ
      4-plex set. Load with \code{data(iTRAQ4)}. }
    \item{\code{iTRAQ5}:}{\code{ReporterIon} object for the iTRAQ
      4-plex set plus the isobaric tag. Load with \code{data(iTRAQ5)}. }
    \item{\code{TMT6}:}{\code{ReporterIon} object for the TMT
      6-plex set. Load with \code{data(TMT6)}. }
    \item{\code{TMT7}:}{\code{ReporterIon} object for the TMT
      6-plex set plus the isobaric tag. Load with \code{data(TMT6)}. }
  }  
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ReporterIons", ...)}.
}


\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"} to identify
      the \code{ReporterIons} instance. }
    \item{\code{reporterNames}:}{Object of class \code{"character"}
      naming each individual reporter of the \code{ReporterIons}
      instance. If not provided explicitely, they are names by
      concatenating the \code{ReporterIons} name and the respective MZ
      values. }
    \item{\code{description}:}{Object of class \code{"character"} to
      describe the \code{ReporterIons} instance. }
    \item{\code{mz}:}{Object of class \code{"numeric"} providing the MZ
      values of the reporter ions. }
    \item{\code{col}:}{Object of class \code{"character"} providing
      colours to highlight the reporters on plots. }
    \item{\code{width}:}{Object of class \code{"numeric"} indicating the
      width around the individual reporter ions MZ values were to search
      for peaks. This is dependent on the mass spectrometer's resolution
      and is used for peak picking when quantifying the reporters. See
      \code{\link{quantify}} for more details about quantification. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}
      indicating the version of the \code{ReporterIons}
      instance. Intended for developer use and debugging.  }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{Versioned}"}, directly.
}

%%\section{Methods}{
%% No methods defined with class "ReporterIons" in the signature.
%%}

\references{
  Ross PL, Huang YN, Marchese JN, Williamson B, Parker K, Hattan S,
  Khainovski N, Pillai S, Dey S, Daniels S, Purkayastha S, Juhasz P,
  Martin S, Bartlet-Jones M, He F, Jacobson A, Pappin DJ.
  "Multiplexed protein quantitation in Saccharomyces cerevisiae using
  amine-reactive isobaric tagging reagents."
  \emph{Mol Cell Proteomics}, 2004 Dec;3(12):1154-69.
  Epub 2004 Sep 22. PubMed PMID: 15385600.
  
  Thompson A, Sch\"{a}fer J, Kuhn K, Kienle S, Schwarz J, Schmidt G,
  Neumann T, Johnstone R, Mohammed AK, Hamon C.
  "Tandem mass tags: a novel quantification strategy for comparative
  analysis of complex protein mixtures by MS/MS."
  \emph{Anal Chem.} 2003 Apr 15;75(8):1895-904. \emph{Erratum} in:
  \emph{Anal Chem.} 2006 Jun 15;78(12):4235. Mohammed, A Karim A [added]
  and  
  \emph{Anal Chem.} 2003 Sep 15;75(18):4942. Johnstone, R [added].
  PubMed PMID: 12713048. 
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
## Code used for the iTRAQ4 set
new("ReporterIons",
    description="4-plex iTRAQ",
    name="iTRAQ4",
    mz=c(114.13,115.13,116.13,117.13),
    col=c("red","green","blue","yellow"),
    width=0.05)
}

\keyword{classes}
