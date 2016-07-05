\name{readMSData}
\alias{readMSData}
\alias{readMSData2}
\concept{MSnExp}

\title{
  Imports mass-spectrometry raw data files as 'MSnExp' instances.
}

\description{
  Reads as set of XML-based mass-spectrometry data files and generates
  an \code{"\linkS4class{MSnExp}"} object. This function uses the
  functionality provided by the \code{mzR} package to access data and
  meta data in \code{mzData}, \code{mzXML} and \code{mzML}.

  \bold{TODO}: document \code{readMSData2}, or whatever it is going to
  be named once stable.
  
}

\usage{
readMSData(files, pdata = NULL, msLevel = 2, verbose = TRUE,
centroided = FALSE, smoothed = FALSE, removePeaks = 0, clean = FALSE,
cache = 1)

readMSData2(files, pdata = NULL, msLevel, verbose = TRUE, centroided =
NA, smoothed = NA)
}

\arguments{
  \item{files}{
    character vector with file names to be read.
  }
  \item{pdata}{
    an object of class \code{"\linkS4class{NAnnotatedDataFrame}"}.
  }
  \item{msLevel}{
    MS level spectra to be read. Use '1' for MS1 spectra or any larger
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
    Numeric indicating caching level. Default is 0 for MS1 and 1 MS2 (or
    higher). Under development.
  }
  \item{verbose}{
   verbosity flag.
  }
}

\details{

  For \code{readMSData} the whole MS data is read from file and kept in
  memory as \code{\linkS4class{Spectrum}} objects within the
  \code{\linkS4class{MSnExp}}'s \code{assayData} slot. To reduce the
  memory footpring especially for large MS 1 data sets it is also
  possible to read only selected information from the MS files and fetch
  the actual spectrum data (i.e. the M/Z and intensity values) only on
  demand from the original data files. This can be achieved with
  \code{readMSData2}. The function returns then an \code{OnDiskMSnExp}
  object instead of a \code{\linkS4class{MSnExp}} object.

}

\value{
  An \code{MSnExp} object for \code{readMSData} and a
  \code{OnDiskMSnExp} object for \code{readMSData2}.
  
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
  \code{"\linkS4class{MSnExp}"} or \code{\link{readMgfData}} to read
  \code{mgf} peak lists.
}

\examples{
file <- dir(system.file(package="MSnbase",dir="extdata"),
            full.name=TRUE,
            pattern="mzXML$")
aa <- readMSData(file)
aa
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ file }
\keyword{ manip }