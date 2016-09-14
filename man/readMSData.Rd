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

}

\usage{

readMSData(files, pdata = NULL, msLevel. = 2, verbose =
isMSnbaseVerbose(), centroided. = NA, smoothed. = NA, cache. = 1)

readMSData2(files, pdata = NULL, msLevel., verbose = isMSnbaseVerbose(),
centroided., smoothed. = NA)

}

\arguments{

  \item{files}{
    A \code{character} vector with file names to be read.
  }

  \item{pdata}{

    An object of class \code{"\linkS4class{NAnnotatedDataFrame}"} or
    \code{NULL}.

  }

  \item{msLevel.}{

    MS level spectra to be read. Use \code{1} for MS1 spectra or any larger
    numeric for MSn spectra. Default is \code{2}. \code{readMSData2}
    supports multiple levels and will, by default, read all the data.

  }

  \item{centroided.}{

    A \code{logical} in \code{readMSData}, indicating whether spectra
    are centroided or not. Default is \code{NA}.

    In \code{readMSData2}, the parameter value is missing and the value
    is set to \code{NA} in the returned object. It can be set for
    different MS levels by a vector of logicals, where the first element
    is for MS 1, the second element is for level 2, ... See
    \code{"\linkS4class{OnDiskMSnExp}"} for an example.

  }

  \item{smoothed.}{

    A \code{logical} indicating whether spectra already smoothed or
    not. Default is \code{NA} in \code{readMSData} and
    \code{readMSData2}.

  }

  \item{cache.}{
    Numeric indicating caching level. Default is 0 for MS1 and 1 MS2 (or
    higher). Under development.
  }
  \item{verbose}{
    Verbosity flag. Default is \code{TRUE}.
  }
}

\details{

  For \code{readMSData} the whole MS data is read from file and kept in
  memory as \code{\linkS4class{Spectrum}} objects within the
  \code{\linkS4class{MSnExp}}'s \code{assayData} slot.

  To reduce the memory footpring especially for large MS 1 data sets it
  is also possible to read only selected information from the MS files
  and fetch the actual spectrum data (i.e. the M/Z and intensity values)
  only on demand from the original data files. This can be achieved with
  \code{readMSData2}. The function returns then an \code{OnDiskMSnExp}
  object instead of a \code{\linkS4class{MSnExp}} object.

}

\value{

  An \code{MSnExp} object for \code{readMSData} and a
  \code{OnDiskMSnExp} object for \code{readMSData2}.

}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}


\seealso{
  \code{"\linkS4class{MSnExp}"} or \code{\link{readMgfData}} to read
  \code{mgf} peak lists.
}

\examples{
file <- dir(system.file(package = "MSnbase", dir = "extdata"),
            full.name = TRUE,
            pattern = "mzXML$")
aa <- readMSData(file)
aa
bb <- readMSData2(file)
bb
}

\keyword{ file }
\keyword{ manip }