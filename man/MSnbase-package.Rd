\name{MSnbase-package}
\alias{MSnbase-package}
\alias{MSnbase}
\docType{package}
\title{MSnbase: Base Functions and Classes for MS-based Proteomics}

\description{
MSnbase provides classes, methods and functions for the plotting,
manipulation and processing of mass spectrometry data.

Important class are \code{"\linkS4class{MSnExp}"} (raw data file),
\code{"\linkS4class{MSnSet}"} (expression data) and
\code{"\linkS4class{ReporterIons}"} (reporter ions for labelled
proteomics).  
 
 Other classes are \code{"\linkS4class{Spectrum}"} and the subclasses
 \code{"\linkS4class{Spectrum1}"} (for MS spectra) and
 \code{"\linkS4class{Spectrum2}"} (for MSMS spectra),
 \code{"\linkS4class{MIAPE} "} (Minimum Information about Proteomics
 Experiments) and MSnProcess (for processing information). 

 Introductory information, use cases and details are available from the
 vignettes, readily accessible with \code{openVignette()}. 

 Complete listing of available documentation with
 \code{library(help="MSnbase")}.
}


\details{
\tabular{ll}{
Package: \tab MSnbase\cr
Version: \tab 0.0.1\cr
Depends: \tab Biobase, methods, ggplot2, plyr\cr
License: \tab Artistic-2.0\cr
Collate: \tab DataClasses.R AllGenerics.R AllMethods.R
formatHandlers.R MSnbase-package.R utils.R
functions-MIAPE.R functions-MSnExp.R functions-MSnProcess.R
functions-MSnQual.R functions-MSnSet.R functions-ReporterIons.R
functions-Spectrum1.R functions-Spectrum2.R functions-Spectrum.R
plotting-MSnExp.R plotting-MSnSet.R plotting-Spectrum1.R
plotting-Spectrum2.R plotting-Spectrum.R\cr
LazyLoad: \tab yes\cr
biocViews: \tab Infrastructure, Bioinformatics, Proteomics\cr
Packaged: \tab 2010-10-08 09:30:44 UTC; lgatto\cr
Built: \tab R 2.13.0; ; 2010-10-08 09:30:45 UTC; unix\cr
}

Index:
\preformatted{
MIAPE-class             The "MIAPE" Class for Storing Proteomics
                        Experiment Information
MSnExp-class            The "MSnExp" Class for MS Data And Meta-Data
MSnProcess-class        The "MSnProcess" Class
MSnSet-class            The "MSnSet" Class for MS Proteomics Expression
                        Data and Meta-Data
ReporterIons-class      The "ReporterIons" Class
Spectrum-class          The "Spectrum" Class
Spectrum1-class         The "Spectrum1" Class for MS1 Spectra
Spectrum2-class         The "Spectrum2" Class for MSMS Spectra
TMT6                    TMT 6-plex set
iTRAQ4                  iTRAQ 4-plex set
}


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
  methods. Several \code{MSnbase} data strucutres are based on
  \code{\link[<pkg>:Biobase-package]{Biobase}} classes.
}


