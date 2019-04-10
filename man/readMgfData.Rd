\name{readMgfData}
\alias{readMgfData}

\title{
  Import mgf files as 'MSnExp' instances.
}

\description{
  Reads a mgf file and generates an \code{"\linkS4class{MSnExp}"}
  object.
}

\usage{
readMgfData(filename, pdata = NULL, centroided = TRUE, smoothed = FALSE,
verbose = isMSnbaseVerbose(), cache = 1)
}


\arguments{
  \item{filename}{
    character vector with file name to be read.
  }
  \item{pdata}{
    an object of class \code{"\linkS4class{AnnotatedDataFrame}"}.
  }
  \item{smoothed}{
    \code{Logical} indicating whether spectra already smoothed or not. Default is
    'FALSE'. Used to initialise \code{"\linkS4class{MSnProcess}"} object
    in \code{processingData} slot.
  }
  \item{centroided}{
    \code{Logical} indicating whether spectra are centroided or not. Default is
    'TRUE'. Used to initialise \code{"\linkS4class{MSnProcess}"} object
    in \code{processingData} slot.
  }
  \item{cache}{
    Numeric indicating caching level. Default is 1. Under development.
  }
  \item{verbose}{
    verbosity flag. 
  }
}

\value{
  An instance of 
}

\author{
  Guangchuang Yu <guangchuangyu@gmail.com> and Laurent Gatto <lg390@cam.ac.uk>
}

\details{
  Note that when reading an mgf file, the original order of the spectra
  is lost. Thus, if the data was originally written to mgf from an
  \code{MSnExp} object using \code{writeMgfData}, although the feature
  names will be identical, the spectra are not as a result of the
  reordering. See example below.
}

\seealso{
  \code{\link{writeMgfData}} method to write the content of
  \code{"\linkS4class{Spectrum}"} or \code{"\linkS4class{MSnExp}"}
  objects to mgf files. Raw data files can also be read with the
  \code{\link{readMSData}} function.
}
  
\examples{
  data(itraqdata)
  writeMgfData(itraqdata, con="itraqdata.mgf", COM="MSnbase itraqdata")
  itraqdata2 <- readMgfData("itraqdata.mgf")
  ## note that the order of the spectra is altered
  ## and precision of some values (precursorMz for instance)
  match(signif(precursorMz(itraqdata2),4),signif(precursorMz(itraqdata),4))
  ## [1]  1 10 11 12 13 14 15 16 17 18 ...
  ## ... but all the precursors are there
  all.equal(sort(precursorMz(itraqdata2)),
            sort(precursorMz(itraqdata)),
            check.attributes=FALSE,
            tolerance=10e-5)
  ## is TRUE
  all.equal(as.data.frame(itraqdata2[[1]]),as.data.frame(itraqdata[[1]]))
  ## is TRUE
  all.equal(as.data.frame(itraqdata2[[3]]),as.data.frame(itraqdata[[11]]))
  ## is TRUE
  f <- dir(system.file(package="MSnbase",dir="extdata"),
           full.name=TRUE,
           pattern="test.mgf")
  (x <- readMgfData(f))
  x[[2]]
  precursorMz(x[[2]])
  precursorIntensity(x[[2]])
  precursorMz(x[[1]])
  precursorIntensity(x[[1]]) ## was not in test.mgf
  scanIndex(x)
}

\keyword{ file }
