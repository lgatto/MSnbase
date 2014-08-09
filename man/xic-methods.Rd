\name{xic-methods}
\docType{methods}
\alias{xic-methods}
\alias{xic,mzRramp-method}
\alias{xic,character-method}
\alias{xic}

\title{ Extracted ion chromatograms }

\description{
  These methods produce an extracted ion chromatogram given
  mass-spectrometry data and an ion to be extracted. In addition to the
  \code{object} (see 'Methods' section below), additional arguments are
  \describe{
    \item{mz}{A \code{numeric} specifying the ion mass to be
      extracted. }
    \item{width}{The M/Z extraction with. Default is 0.5.}
    \item{rtlim}{The retention time limit to be displayed. When missing
      (default), the complete range of the matching extracted ions is
      plotted. }
    \item{npeaks}{The number of peaks to be annotated. Default is 3.}  
    \item{charge}{The charge of the ion to be extracted. This value is
      optional; when provided, the mass of the ion (\code{mz} above)
      will first be divided by its \code{charge} before extraction. }
    \item{clean}{A \code{logical}, defining of the XIC data should be
      cleaned before plotting. Default if \code{TRUE}. See
      \code{\link{clean}} for details. } 
    \item{legend}{A \code{logical} defining if the figure should be
      annotated. }
    \item{plot}{A \code{logical} defining if the plot should be
      rendered. }
    \item{points}{A \code{logical} specifying if points should be added
      to for each individual MS spectrum. Annotated peaks are coloured in
      orange and MS2 spectra with precursor M/Z matching the extracted
      ion are highlighted in red. }
    \item{hd}{The header \code{data.frame} corresponding to the
      \code{object} input. \code{hd} is generated automatically and can
      generally be omitted. }
    \item{...}{Additional arguments passed to the \code{plot} function.}
  }
  Note: the methods are currently not vectorised.

  xcms::plotEIC provides a similar functionality.
}

\value{
  The methods invisibly return the \code{data.frame} with the the XIC
  intensities (column \code{int}), retention times (column \code{rt})
  and extracted M/Z values (column \code{mz}) used to generate the
  figure.  
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "character")}}{ Plots the XIC for the
      mass-spectrometry data stored in the \code{object} file. The file
      format must be support by \code{mzR}. See \code{mzR::openMSfile} for
      details. 
    }
    \item{\code{signature(object = "mzRramp")}}{ Plots the XIC
      for the \code{mzRramp} instance. See the \code{mzR} package for
      details.   
    }
  }
}

\examples{
\dontrun{
    library("RforProteomics")
    (f <- getPXD000001mzXML())
    ms <- openMSfile(f)
    x <- xic(ms, mz = 636.925, width = 0.01)
    head(x)
    xic(ms, mz = 636.925, width = 0.01,
        npeaks = 1,
        rtlim = c(2000, 2300))
}
}

\keyword{methods}

