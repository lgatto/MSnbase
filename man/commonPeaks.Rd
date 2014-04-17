\name{commonPeaks}
\alias{commonPeaks}
\title{common peaks}
\usage{
commonPeaks(x, y, method = c("highest", "closest"), tolerance = 2.5e-05,
  relative = TRUE)
}
\arguments{
  \item{x}{spectrum1 (MSnbase::Spectrum), to be matched}

  \item{y}{spectrum2 (MSnbase::Spectrum), to match against}

  \item{method}{for duplicated matches use highest/closest
  intensity/mz}

  \item{tolerance}{double, allowed deviation}

  \item{relative}{relative (or absolute) deviation}
}
\value{
logical, TRUE for common peaks in y
}
\description{
common peaks
}

