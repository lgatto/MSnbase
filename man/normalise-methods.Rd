\name{normalise-methods}
\docType{methods}
\alias{normalise-methods}
\alias{normalise,Spectrum-method}
\alias{normalise,Spectrum2-method}
\alias{normalise,MSnSet-method}
\alias{normalise,MSnExp-method}
\alias{normalize-methods}
\alias{normalize,Spectrum-method}
\alias{normalize,Spectrum2-method}
\alias{normalize,MSnSet-method}
\alias{normalize,MSnExp-method}
\alias{scale,MSnSet-method}
\alias{normalise}
\alias{normalize}


\title{ Normalisation of \code{MSnExp}, \code{MSnSet} and
  \code{Spectrum} objects } 

\description{
  The \code{normalise} method (also available as \code{normalize})
  performs basic normalisation on spectra 
  intensities of single spectra (\code{"\linkS4class{Spectrum}"} or 
  \code{"\linkS4class{Spectrum2}"} objects),
  whole experiments (\code{"\linkS4class{MSnExp}"} objects) or
  quantified expression data (\code{"\linkS4class{MSnSet}"} objects).

  Raw spectra and experiments are normalised using \code{max} or
  \code{sum} only. For MSMS spectra could be normalised to their 
  \code{precursor} additionally. Each peak intensity is divided by the
  highest intensity in the spectrum, the sum of intensities or the intensity 
  of the precursor.
  These methods aim at facilitating relative peaks heights between
  different spectra.

  The \code{method} parameter for \code{"\linkS4class{MSnSet}"} can be
  one of \code{sum}, \code{max}, \code{quantiles}, \code{center.mean},
  \code{center.median}, \code{.median}, \code{quantiles.robust} or
  \code{vsn}.  For \code{sum} and \code{max}, each feature's reporter
  intensity is divided by the maximum or the sum respectively. These two
  methods are applied along the features (rows).
  
  \code{center.mean} and \code{center.median} translate the respective
  sample (column) intensities according to the column mean or
  median. \code{diff.median} translates all samples (columns) so that
  they all match the grand median. Using \code{quantiles} or
  \code{quantiles.robust} applies (robust) quantile normalisation, as
  implemented in \code{normalize.quantiles} and
  \code{normalize.quantiles.robust} of the \code{preprocessCore}
  package.  \code{vsn} uses the \code{\link{vsn2}} function from the
  \code{vsn} package.  Note that the latter also glog-transforms the
  intensities.  See respective manuals for more details and function
  arguments.
  
  A \code{scale} method, mimicking the base \code{scale} method exists
  for \code{"\linkS4class{MSnSet}"} instances. See
  \code{?base::\link[base]{scale}} for details. 
}

\arguments{
  \item{object}{ An object of class \code{"\linkS4class{Spectrum}"}, 
    \code{"\linkS4class{Spectrum2}"}, \code{"\linkS4class{MSnExp}"} or
    \code{"\linkS4class{MSnSet}"}.
  }
  \item{method}{
    A character vector of length one that describes how to normalise
    the object. See description for details.
  }
  \item{...}{ Additional arguments passed to the normalisation function.
  }
}

\section{Methods}{
  The \code{normalise} methods:
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
    \item{\code{signature(object = "Spectrum2", method = "character",
        precursorIntensity)}}{
      Normalises the \code{object} peak intensities using
      \code{method}. If \code{method == "precursor"}, 
      \code{precursorIntensity} allows to specify the intensity of the
      precursor manually. }
  }
  The \code{scale} method:
  \describe{
    \item{\code{signature(x = "MSnSet", center = "logical", scale =
	"logical")}}{ See \code{?base::\link[base]{scale}}. }
  }
}

\examples{
## quantifying full experiment
data(msnset)
msnset.nrm <- normalise(msnset, "quantiles")
msnset.nrm
}


\keyword{methods}
