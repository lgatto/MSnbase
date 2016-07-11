\name{bin-methods}
\docType{methods}
\alias{bin-methods}
% \alias{bin,MSnExp-method}
% \alias{bin,Spectrum-method}
\alias{bin}


\title{ Bin 'MSnExp' or 'Spectrum' instances }

\description{
  This method aggregates individual spectra (\code{Spectrum} instances)
  or whole experiments (\code{MSnExp} instances) into discrete bins. All
  intensity values which belong to the same bin are summed together.
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", binSize = "numeric", verbose =
	"logical")}}{ Bins all spectra in an \code{MSnExp}
      object. Use \code{binSize} to control the size of a bin
      (in Dalton, default is \code{1}).
      Displays a control bar if verbose set to \code{TRUE}
      (default). Returns a binned \code{MSnExp} instance.
    }

    \item{\code{signature(object = "Spectrum", binSize = "numeric",
        breaks = "numeric", msLevel. = "numeric")}}{ Bin the
        \code{Spectrum} object.  Use \code{binSize} to control the size
        of a bin (in Dalton, default is \code{1}).  Similar to
        \code{\link[graphics]{hist}} you could use \code{breaks} to
        specify the breakpoints between m/z bins. \code{msLevel.}
        defines the level of the spectrum, and if \code{msLevel(object)
        != msLevel.}, cleaning is ignored. Only relevant when called
        from \code{OnDiskMSnExp} and is only relevant for developers.
      
	Returns a binned \code{Spectrum} instance.
    }
  }
}

\author{
  Sebastian Gibb <mail@sebastiangibb.de>
}

\seealso{
  \code{\link{clean}}, \code{\link{pickPeaks}}, \code{\link{smooth}},
  \code{\link{removePeaks}} and \code{\link{trimMz}}
  for other spectra processing methods.
}

\examples{
s <- new("Spectrum2", mz=1:10, intensity=1:10)
intensity(s)
intensity(bin(s, binSize=2))

data(itraqdata)
sum(peaksCount(itraqdata))
itraqdata2 <- bin(itraqdata, binSize=2)
sum(peaksCount(itraqdata2))
processingData(itraqdata2)
}


\keyword{methods}
