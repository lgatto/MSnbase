\name{readMzXMLData}
\alias{readMzXMLData}
\concept{MSnExp}

\title{
  Import mzXML files as 'MSnExp' instances.
}

\description{
  Reads mzXML files and generates an \code{"\linkS4class{MSnExp}"}
  object.

  Note: the preferred function to import raw data to generate
  \code{MSnExp} objects is \code{\link{readMSData}}. 
}

\usage{
readMzXMLData(files, pdata = NULL, msLevel = 2, verbose = TRUE,
centroided = FALSE, smoothed = FALSE, removePeaks = 0, clean = FALSE,
cache = 1) 
}

\arguments{
  \item{files}{
    character vector with file names to be read.
  }
  \item{pdata}{
    an object of class \code{"\linkS4class{NAnnotatedDataFrame}"}.
  }
  \item{msLevel}{
    MS level spectra to be read. Use '1' for MS1 spectra of any larger
    numeric for MSn spectra. Default is '2'.
  }
  \item{centroided}{
    \code{Logical} indicating whether spectra are centroided or not. Default is
    'FALSE'. Used to initialise \code{"\linkS4class{MSnProcess}"} object
    in \code{processingData} slot.
  }
  \item{smoothed}{
    \code{Logical} indicating whether spectra already smoothed or not. Default is
    'FALSE'. Used to initialise \code{"\linkS4class{MSnProcess}"} object
    in \code{processingData} slot.
  }
  \item{removePeaks}{
    If > 0 (default), all peaks less or equal then value will set to
    0. See \code{\link{removePeaks}} for more details and examples.
  }
  \item{clean}{
    Logical indicating whether 0 intensity peaks should be discarded
    from spectra. Useful is removePeaks is set. Default is 'FALSE'. See
    \code{\link{clean}} for more details and examples. 
  }
  \item{cache}{
    Numeric indicating caching level. Default is 1. Under development.   
  }
  \item{verbose}{
   verbosity flag. 
  }
}

% \details{
%   %%  ~~ If necessary, more details than the description above ~~
% }

\value{
  An \code{"\linkS4class{MSnExp}"} object.
}

% \references{
%   %% ~put references to the literature/web site here ~
% }

\author{
Laurent Gatto <lg390@cam.ac.uk>
}

% \note{
% %%  ~~further notes~~
% }



\seealso{
  \code{"\linkS4class{MSnExp}"} as well as the \code{\link{readMSData}}
  function to read more formats.
}

\examples{
file <- dir(system.file(package="MSnbase",dir="extdata"),
            full.name=TRUE,
            pattern="mzXML$")
aa <- readMzXMLData(file)
aa
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ manip }