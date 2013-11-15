\name{readMzTabData}
\alias{readMzTabData}
\title{Read an 'mzTab' file}
\usage{
  readMzTabData(file, what = c("PRT", "PEP"),
    verbose = TRUE)
}
\arguments{
  \item{file}{A \code{character} with the \code{mzTab} file
  to be read in.}

  \item{what}{One of \code{"PRT"} or \code{"PEP"}, defining
  which of protein of peptide section should be parse. The
  metadata section, when available, is always used to
  populate the \code{experimentData} slot.}

  \item{verbose}{Produce verbose output.}
}
\value{
  An instance of class \code{MSnSet}.
}
\description{
  This function can be used to create a
  \code{"\linkS4class{MSnSet}"} by reading and parsing an
  \code{mzTab} file. The metadata section is always used to
  populate the \code{MSnSet}'s \code{experimentData} slot.
}
\examples{
testfile <- "http://mztab.googlecode.com/svn/legacy/jmztab-1.0/examples/mztab_itraq_example.txt"
prot <- readMzTabData(testfile, "PRT")
prot
pep <- readMzTabData(testfile, "PEP")
pep
}
\author{
  Laurent Gatto
}
\seealso{
  \code{\link{writeMzTabData}} to save an
  \code{"\linkS4class{MSnSet}"} as an \code{mzTab} file.
}

