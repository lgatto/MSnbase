\name{writeMgfData-methods}
\docType{methods}
\alias{writeMgfData}
\alias{writeMgfData-methods}
\alias{writeMgfData,MSnExp-method}
\alias{writeMgfData,Spectrum-method}

\title{ Write an experiment or spectrum to an mgf file }

\description{
  Methods \code{writeMgfData} write individual
  \code{"\linkS4class{Spectrum}"} instances of whole
  \code{"\linkS4class{MSnExp}"} experiments to a file
  in Mascot Generic Format (mgf) (see 
  \url{http://www.matrixscience.com/help/data_file_help.html}
  for more details). Function \code{readMgfData} read spectra from and
  mgf file and creates an \code{"\linkS4class{MSnExp}"} object.
  
}

\arguments{
  \item{object}{An instance of class \code{"\linkS4class{Spectrum}"} or
    \code{"\linkS4class{MSnExp}"}. }
  \item{con}{A valid \code{connection} or a \code{character} string with
    the name of the file to save the object. In case of the latter, a
    \code{file} connection is created. If not specified, 'spectrum.mgf' 
    or 'experiment.mgf' are used depending on the class of
    \code{object}. Note that existing files are overwritted. }
  \item{COM}{Optional character vector with the value for the 'COM' field. }
  \item{TITLE}{Optional character vector with the value for the spectrum
    'TITLE' field. Not applicable for experiments.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "MSnExp")}}{Writes the full
      exeriment to an mgf file. }
    \item{\code{signature(object = "Spectrum")}}{Writes an individual
      spectrum to an mgf file. }
  }
}

\seealso{
  \code{\link{readMgfData}} function to read data from and mgf file.
}

\details{
  Note that when reading an mgf file, the original order of the spectra
  is lost. Thus, if the data was originally written to mgf from an
  \code{MSnExp} object using \code{writeMgfData}, although the feature
  names will be identical, the spectra are not as a result of the
  reordering. See example below.
}

\examples{
data(itraqdata)

f <- tempfile()

writeMgfData(itraqdata, con = f)

itraqdata2 <- readMgfData(f)

## note that the order of the spectra and precision of some values
## (precursorMz for instance) are altered
match(signif(precursorMz(itraqdata2),4),
      signif(precursorMz(itraqdata),4))

## [1]  1 10 11 12 13 14 15 16 17 18 ...
## ... but all the precursors are there
all.equal(sort(precursorMz(itraqdata2)),
          sort(precursorMz(itraqdata)),
          check.attributes = FALSE,
          tolerance = 10e-5)

all.equal(as.data.frame(itraqdata2[[1]]),
          as.data.frame(itraqdata[[1]]))

all.equal(as.data.frame(itraqdata2[[3]]),
          as.data.frame(itraqdata[[11]]))

all(featureNames(itraqdata2) == featureNames(itraqdata))
}


\keyword{file}
\keyword{methods}
