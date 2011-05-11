\name{writeMgfData-methods}
\docType{methods}
\alias{writeMgfData}
\alias{writeMgfData-methods}
\alias{writeMgfData,MSnExp-method}
\alias{writeMgfData,Spectrum-method}

\title{ Writes an experiment or spectrum to an mgf file }

\description{
  Methods \code{writeMgfData} write individual
  \code{"\linkS4class{Spectrum}"} instances of whole
  \code{"\linkS4class{MSnExp}"} experiments to a file
  in Mascot Generic Format (mgf) (see 
  \url{http://www.matrixscience.com/help/data_file_help.html}
  for more details)
  
}

\arguments{
  \item{object}{An instance of class \code{"\linkS4class{Spectrum}"} or
    \code{"\linkS4class{MSnExp}"}. }
  \item{filename}{A character string with the name of the file to save
    the object. If not specified, 'spectrum.mgf' or 'experiment.mgf' are
    used depending on the class of \code{object}. }
  \item{COM}{Optional character vector with the value for the 'COM' field. }
  \item{TITLE}{Optional character vector with the value for the spectrrm
    'TITLE' field. Not applicable for experiments.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp")}}{Writes the full
      exeriment to an mgf file. }
    \item{\code{signature(object = "Spectrum")}}{Writes an individual
      spectrum to an mgf file. }
  }
}

\keyword{methods}
