\name{formatRt}
\alias{formatRt}


\title{
  Format Retention Time
}

\description{
  Converts seconds to min:sec format
}

\usage{
formatRt(rt)
}

\arguments{
  \item{rt}{ retention in in seconds }
}

\details{
  This function is used to convert retention times, expressed in seconds,
  in the more human friendly format mm:sec.
}

\value{
  a character string in mm:ss format
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
formatRt(1524)
}

\keyword{ chron }
\keyword{ utilities }
