\name{MSnbase-package}
\alias{MSnbase-package}
\alias{MSnbase}
\docType{package}
\title{MSnbase: Base Functions and Classes for MS-based Proteomics}

\description{

  MSnbase provides classes, methods and functions for visualisation,
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

  If you have questions, want to rebort a bug or share suggestions,
  please file an issue at
  \url{https://github.com/lgatto/MSnbase/issues}, contact me directly or
  ask a question on the Bioconductor support forum
  \url{https://support.bioconductor.org/}.
  
}


\author{
  Laurent Gatto
  
  Maintainer: Laurent Gatto <lg390@cam.ac.uk>

  See the DESCRIPTION file for a complete list of contributors.
}

\references{
  Laurent Gatto and Kathryn S. Lilley,
  MSnbase - an R/Bioconductor package for isobaric tagged mass
  spectrometry data visualization, processing and quantitation,
  Bioinformatics 28(2), 288-289 (2012).
  
  Gatto L. and Lilley K.S., Towards reproducible MSMS data
  preprocessing, quality control and quantification. BSPR/EBI Proteomics
  Meeting,  Hinxton, United Kingdom, 13-15 July 2010,
  \url{http://dx.doi.org/10.1038/npre.2010.5010.1}. 
}

\keyword{ package }

\seealso{
  Introductory information, use cases and details are available from the
  vignettes:
  \itemize{
    \item The demo vignette describe an use-case using a dummy data set
    provided with the package. It can be accessed with
    \code{vignette("MSnbase-demo", package = "MSnbase")}. 

    \item The development vignette describes the classes implemented in
    \code{MSnbase} and can be accessed with
    \code{vignette("MSnbase-development", package = "MSnbase")}.

    \item Details about input/outcupt capabilities and formats
    can be found in 
    \code{vignette("MSnbase-io", package = "MSnbase")}.
    
  }
  
  Complete listing of available documentation with
  \code{library(help = "MSnbase")}.
}


