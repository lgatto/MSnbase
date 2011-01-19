\name{MSnbase-package}
\alias{MSnbase-package}
\alias{MSnbase}
\docType{package}
\title{MSnbase: Base Functions and Classes for MS-based Proteomics}

\description{
MSnbase provides classes, methods and functions for the plotting,
manipulation and processing of mass spectrometry data.

Important class are \code{"\linkS4class{MSnExp}"} (raw data file),
\code{"\linkS4class{MSnSet}"} (quantitation data) and
\code{"\linkS4class{ReporterIons}"} (reporter ions for labelled
proteomics).  

Other classes are \code{"\linkS4class{Spectrum}"} and the subclasses
\code{"\linkS4class{Spectrum1}"} (for MS spectra) and
\code{"\linkS4class{Spectrum2}"} (for MSMS spectra),
\code{"\linkS4class{MIAPE}"} (Minimum Information about Proteomics
Experiments) and \code{"\linkS4class{MSnProcess}"} (for processing
information). These should however not be of direct utility to users. 

Introductory information, use cases and details are available from the
vignettes, readily accessible with \code{openVignette()}. 

Complete listing of available documentation with
\code{library(help="MSnbase")}.
}


\details{
  see 'DESCRIPTION' file.
}


\author{
  Laurent Gatto
  Maintainer: Laurent Gatto <lg390@cam.ac.uk>
}

\references{
  Package sources: \url{http://github.com/lgatto/MSnbase}
}

\keyword{ package }

\seealso{
  See also the package vignette for use-cases and details about the
  methods. Several \code{MSnbase} data structures are based on
  \code{\link[<pkg>:Biobase-package]{Biobase}} classes.
}


