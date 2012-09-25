\name{extractPrecSpectra-methods}
\docType{methods}
\alias{extractPrecSpectra-methods}
\alias{extractPrecSpectra}


\title{ Extracts precursor-specific spectra from an 'MSnExp' object }

\description{
  Extracts the MSMS spectra that originate from the precursor(s) having
  the same MZ value as defined in the\code{prec} argument.

  A warning will be issued of one or several of the precursor MZ values
  in \code{prec} are absent in the experiment precursor MZ values (i.e
  in \code{precursorMz(object)}).
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp", prec = "numeric")}}{
      Returns an \code{"\linkS4class{MSnExp}"} containing MSMS spectra
      whose precursor MZ values are in \code{prec}. 
    }
  }
}

\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

% \seealso{
%   \code{\link{extractSpectra}}
% }

\examples{
file <- dir(system.file(package="MSnbase",dir="extdata"),
            full.name=TRUE,pattern="mzXML$")
aa <- readMSData(file,verbose=FALSE)
my.prec <- precursorMz(aa)[1]
my.prec
bb <- extractPrecSpectra(aa,my.prec)
precursorMz(bb)
processingData(bb)
}


\keyword{methods}

