\name{matchPeaks}
\alias{matchPeaks}
\title{similar to base::match but with tolerance
if there are duplicated matches the highest/closest is choosen}
\usage{
matchPeaks(x, y, method = c("highest", "closest"), tolerance = 2.5e-05,
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
integer vector of the same length as "x" representing the
position in "y"
}
\description{
similar to base::match but with tolerance if there are
duplicated matches the highest/closest is choosen
}

