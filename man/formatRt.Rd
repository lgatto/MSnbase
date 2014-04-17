\name{formatRt}
\alias{formatRt}


\title{
  Format Retention Time
}

\description{
  Converts seconds to/from 'min:sec' format
}

\usage{
formatRt(rt)
}

\arguments{
  \item{rt}{ retention in seconds (\code{numeric}) or "mm:sec"
  (\code{character})}.
}


\details{
  This function is used to convert retention times. Conversion is
  seconds to/from the more human friendly format "mm:sec".
}

\value{
  A vector of same length as \code{rt}.
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
formatRt(1524)
formatRt("25:24")
}

\keyword{ chron }
\keyword{ utilities }
