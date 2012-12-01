\name{normalise-methods}
\docType{methods}
\alias{normalise-methods}
\alias{normalise,Spectrum-method}
\alias{normalise,MSnSet-method}
\alias{normalise,MSnExp-method}
\alias{normalize-methods}
\alias{normalize,Spectrum-method}
\alias{normalize,MSnSet-method}
\alias{normalize,MSnExp-method}
\alias{normalise}
\alias{normalize}


\title{ Normalisation of \code{MSnExp}, \code{MSnSet} and
  \code{Spectrum} objects } 

\description{
  The \code{normalise} method performs basic normalisation on spectra
  intensities of single spectra (\code{"\linkS4class{Spectrum}"}
  objects), whole experiments (\code{"\linkS4class{MSnExp}"} objects) or
  quantified expression data (\code{"\linkS4class{MSnSet}"} objects).

  Raw spectra and experiments are normalised using \code{max} or
  \code{sum} only. Each peak intensity is divided by the highest
  intensity in the spectrum or the sum of intensities.
  These methods aim at facilitating relative peaks heights between
  different spectra.

  The \code{method} parameter for \code{"\linkS4class{MSnSet}"} can be
  one of \code{sum}, \code{max}, \code{quantiles},
  \code{quantiles.robust} or \code{vsn}.
  For \code{sum} and \code{max}, each feature's reporter intensity is
  divided by the maximum of the sum respectively. These two methods are 
  applied along the features (rows).
  \code{scale.mean} and \code{scale.median} scale the respecitve sample
  (column) intensities according to their mean or median.  
  Using \code{quantiles} or  \code{quantiles.robust} uses (robust) 
  quantile normalisation, as implemented in
  \code{\link{normalize.quantiles}} and
  \code{\link{normalize.quantiles.robust}} of the \code{preprocessCore} 
  package.
  \code{vsn} uses the \code{\link{vsn2}} from the \code{vsn} package.
  Note that the latter also glog-transforms the intensities.
  See respective manuals for more details and function arguments. 
}

\arguments{
  \item{object}{ An object of class \code{"\linkS4class{Spectrum}"}, 
    \code{"\linkS4class{MSnExp}"} or \code{"\linkS4class{MSnSet}"}.
  }
  \item{method}{
    A character vector of length one that describes how to normalise
    the object. See description for details.
  }
  \item{...}{ Additional arguments passed to the normalisation function.
  }
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnSet", method = "character")}}{
      Normalises the \code{object} reporter ions intensities using
      \code{method}. }
    \item{\code{signature(object = "MSnExp", method = "character")}}{
      Normalises the \code{object} peak intensities using
      \code{method}. }
    \item{\code{signature(object = "Spectrum", method = "character")}}{
      Normalises the \code{object} peak intensities using
      \code{method}. }
  }
}

\examples{
## quantifying full experiment
data(itraqdata)
qnt <- quantify(itraqdata,method="trap",reporters=iTRAQ4)
qnt.nrm <- normalise(qnt,"quantiles")
qnt.nrm
}


\keyword{methods}
