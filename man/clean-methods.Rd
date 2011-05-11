\name{clean-methods}
\docType{methods}
\alias{clean-methods}
% \alias{clean,MSnExp-method}
% \alias{clean,Spectrum-method}
\alias{clean}


\title{ Cleans 'MSnExp' or 'Spectrum' instances }

\description{
  This method cleans out individual spectra (\code{Spectrum} instances)
  or whole experiments (\code{MSnExp} instances) of 0-intensity
  peaks. Origninal 0-intensity values are retained only around peaks. If
  more than two 0's were separating two peaks, only the first and last
  ones, those directly adjacent to the peak ranges are kept. If two
  peaks are separated by only one 0-intensity value, it is retained. An
  illustrative example is shown below. 
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", verbose = "logical")}}{
      Cleans all spectra in \code{MSnExp} object. Displays a control bar
      if verbose set to \code{TRUE} (default). Returns a cleaned
      \code{MSnExp} instance.
    }
    \item{\code{signature(object = "Spectrum")}}{
      Cleans the \code{Spectrum} object. Returns a cleaned
      \code{Spectrum} instance.
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
sp2 <- clean(sp1)
intensity(sp1)
intensity(sp2)

file <- dir(system.file(package="MSnbase",dir="extdata"),
            full.name=TRUE,pattern="mzXML$")
aa <- readMzXMLData(file,verbose=FALSE)
bb <- clean(aa,verbose=FALSE)
sum(peaksCount(aa))
sum(peaksCount(bb))
processingData(bb)
}


\keyword{methods}
