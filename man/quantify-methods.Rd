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
    the area under the peak(s), the maximum of the peak(s) or the sum of
    all intensities of the peak(s). 
  }
  \item{reporters}{
    An object of class \code{"\linkS4class{ReporterIons}"} that defines
    the peak(s) to be quantified.
  }
  \item{strict}{
    If strict is 'FALSE' (default), the quantitation is performed
    using data points along the entire width of a peak. If strict is set
    to 'TRUE', once the apex(es) is/are identified, only data points
    within apex +/- width of reporter (see 
    \code{"\linkS4class{ReporterIons}"}) are used for quantitation.
  }
  \item{verbose}{
    Verbose of the output (only for \code{MSnExp} objects).
  }
}

\details{
  \code{"\linkS4class{ReporterIons}"} define specific MZ at
  which peaks are expected and a window around that MZ value.
  A peak of interest is searched for in that window.
  Since version 1.1.2, warnings are not thrown anymore in case no data
  is found in that region or if the peak extends outside the window.
  This can be checked manually after quantitation,
  by inspecting the quantitation data (using the \code{exprs} accessor)
  for \code{NA} values or by comaring the \code{lowerMz} and
  \code{upperMz} columns in the \code{"\linkS4class{MSnSet}"}
  \code{qual} slot against the respective expected \code{mz(reporters)}
  +/- \code{width(reporters)}.
  This is illustrated in the example below.
  
  Once the range of the curve is found, quantification is performed.
  If no data points are found in the expected region, 'NA' is   
  returned for the reporter peak MZ.
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
## quantifying full experiment
file <- dir(system.file(package="MSnbase",dir="extdata"),full.name=TRUE,pattern="mzXML$")
aa <- readMzXMLData(file,verbose=FALSE)
msnset <- quantify(aa,method="trap",reporters=iTRAQ4)
msnset
## checking for non-quantified peaks
sum(is.na(exprs(msnset)))
## checking for peaks outside of first reporter width
lowerMz.114 <- qual(msnset)$lowerMz[qual(msnset)$reporter=="114.1"]
upperMz.114 <- qual(msnset)$upperMz[qual(msnset)$reporter=="114.1"]
any(lowerMz.114 < mz(iTRAQ4[1]) - width(iTRAQ4[1]))
any(upperMz.114 > mz(iTRAQ4[1]) + width(iTRAQ4[1]))
## quantifying single spectrum
qty <- quantify(aa[[1]],method="trap",iTRAQ4[1])
qty$peakQuant
qty$curveStats
}


\keyword{methods}

