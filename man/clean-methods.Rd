\name{clean-methods}
\docType{methods}
\alias{clean-methods}
% \alias{clean,MSnExp-method}
% \alias{clean,Spectrum-method}
\alias{clean}


\title{ Clean 'MSnExp', 'Spectrum' or 'Chromatogram' instances }

\description{
  This method cleans out individual spectra (\code{Spectrum} instances),
  chromatograms (\code{\linkS4class{Chromatogram}} instances)
  or whole experiments (\code{MSnExp} instances) of 0-intensity
  peaks. Unless \code{all} is set to \code{FALSE}, original 0-intensity
  values are retained only around peaks. If more than two 0's were
  separating two peaks, only the first and last ones, those directly
  adjacent to the peak ranges are kept. If two peaks are separated by
  only one 0-intensity value, it is retained. An illustrative example is
  shown below.  
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", all = "logical", verbose =
	"logical")}}{ Cleans all spectra in \code{MSnExp}
      object. Displays a control bar if verbose set to \code{TRUE}
      (default). Returns a cleaned \code{MSnExp} instance.
    }

    \item{\code{signature(object = "Spectrum", all = "logical",
	msLevel. = "numeric")}}{ Cleans the \code{Spectrum}
      object. Returns a cleaned \code{Spectrum} instance. If \code{all
	= TRUE}, then all zeros are removed. \code{msLevel.} defines the
      level of the spectrum, and if \code{msLevel(object) !=
	msLevel.}, cleaning is ignored. Only relevant when called from
      \code{OnDiskMSnExp} and is only relevant for developers.
    }

    \item{\code{signature(object = "Chromatogram", all = "logical",
	na.rm = "logical")}}{
      Cleans the \code{\linkS4class{Chromatogram}} instance and returns
      a cleaned \code{\linkS4class{Chromatogram}} object. If
      \code{na.rm} is \code{TRUE} (default is \code{FALSE}) all
      \code{NA} intensities are removed before cleaning the chromatogram.
    }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\seealso{
  \code{\link{removePeaks}} and \code{\link{trimMz}} for other spectra
  processing methods.
}

\examples{
int <- c(1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0)
sp1 <- new("Spectrum2",
           intensity=int,
           mz=1:length(int))
sp2 <- clean(sp1) ## default is all=FALSE
intensity(sp1)
intensity(sp2)
intensity(clean(sp1, all = TRUE))

mz(sp1)
mz(sp2)
mz(clean(sp1, all = TRUE))

data(itraqdata)
itraqdata2 <- clean(itraqdata)
sum(peaksCount(itraqdata))
sum(peaksCount(itraqdata2))
processingData(itraqdata2)

## Create a simple Chromatogram object
chr <- Chromatogram(rtime = 1:12,
                    intensity = c(0, 0, 20, 0, 0, 0, 123, 124343, 3432, 0, 0, 0))

## Remove 0-intensity values keeping those adjacent to peaks
chr <- clean(chr)
intensity(chr)

## Remove all 0-intensity values
chr <- clean(chr, all = TRUE)
intensity(chr)

## Clean a Chromatogram with NAs.
chr <- Chromatogram(rtime = 1:12,
                    intensity = c(0, 0, 20, NA, NA, 0, 123, 124343, 3432, 0, 0, 0))
chr <- clean(chr, all = FALSE, na.rm = TRUE)
intensity(chr)

}


\keyword{methods}
