\name{relaxedMatch}
\alias{relaxedMatch}
\title{similar to base::match but with tolerance in contrast to match the
'table' argument must be sorted and it is possible to get
duplicated matches}
\usage{
relaxedMatch(x, table, nomatch = NA_integer_, tolerance = 2.5e-05,
  relative = TRUE)
}
\arguments{
  \item{x}{values to be matched}

  \item{table}{the values to be matched against, must be
  sorted!}

  \item{nomatch}{the value to be returned in the case when
  no match is found}

  \item{tolerance}{double, allowed deviation}

  \item{relative}{relative (or absolute) deviation}
}
\value{
numeric vector of the same length as "x"
}
\description{
similar to base::match but with tolerance in contrast to
match the 'table' argument must be sorted and it is
possible to get duplicated matches
}

