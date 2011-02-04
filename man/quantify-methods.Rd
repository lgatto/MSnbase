\name{quantify-methods}
\docType{methods}
\alias{quantify-methods}
\alias{quantify}

\title{ Quantifies 'MSnExp' and 'Spectrum' objects }


\description{
  This method quantifies individual \code{"\linkS4class{Spectrum}"}
  objects or full \code{"\linkS4class{MSnExp}"} experiments. Current
  implementation quantify specific peaks in individual spectra. The
  peaks of interest are defined by the \code{reporters} parameter using
  one of the possible \code{method} methods. 
  This essentially qualifies for MSMS quantitation using isobaric tags,
  altough any peak can be defined in \code{reporters}. 
}

\arguments{
  \item{object}{ Objects of class \code{"\linkS4class{Spectrum}"} or
    \code{"\linkS4class{MSnExp}"}. 
  }
  \item{method}{
    Peak quantitation method. One of, possibly abreviated
    "trapezoidation", "max", or "sum". These methods return respectively
    the area under the peak, the maximum of the peak or the sum of all
    intensities of the peaks. 
  }
  \item{reporters}{
    An object of class \code{"\linkS4class{ReporterIons}"} that defines
    the peaks to be quantified.
  }
  \item{verbose}{
    Verbose of the output (only for \code{MSnExp} objects).
  }
}

\details{
  Peaks \code{"\linkS4class{ReporterIons}"} defined a specific MZ at
  which the peak is expected and a window around that MZ value. A peak
  is searched for in that window and warnings are in case no peak data
  is found in that region or the reporter peak extends outside the
  window. Once the range of the curve is found, quantification is
  performed. If no data points are found in the expected region, 'NA' is
  returned at the expected reporter peak MZ.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", method = "character",
	reporters = "ReporterIons", verbose = "logical")}}{
      Quantifies peaks defined in \code{reporters} using \code{method}
      in all spectra of the \code{MSnExp} object. If verbose is set to
      \code{TRUE}, a progress bar will be displayed.
      
      An object of class \code{"\linkS4class{MSnSet}"} is returned
      containing the quantified feature expression and all meta data
      inherited from the \code{MSnExp} argument. 
    }
    \item{\code{signature(object = "Spectrum", method = "character",
	reporters = "ReporterIons")}}{
      Quantifies peaks defined in \code{reporters} using \code{method}
      in the \code{Spectrum} object.

      A list of length 2 will be returned. The first element, named
      \code{peakQuant}, is a 'numeric' of length equal to
      \code{length(reporters)} with quantitation of the reporter peaks
      using \code{method}.

      The second element, names \code{curveStats}, is a 'data.frame'
      of dimension \code{length(reporters)} times 7 giving, for each
      reporter curve parameters: maximum intensity ('maxInt'), number
      of maxima ('nMaxInt'), number of data points defined the curve
      ('baseLength'), lower and upper MZ values for the curve
      ('lowerMz' and 'upperMz'), reporter ('reporter') and precursor
      MZ value ('precursor') when available. 	
    }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
aa <- readMzXMLData(file,verbose=FALSE)
data(iTRAQ4)
msnset <- quantify(aa,method="trap",reporters=iTRAQ4)
msnset
qty64 <- quantify(aa[["X64"]],method="trap",iTRAQ4[1])
qty64$peakQuant
qty64$curveStats
}


\keyword{methods}

