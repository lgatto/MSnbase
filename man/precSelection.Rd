\name{precSelection}
\alias{precSelection}
\alias{precSelectionTable}


\title{ Number of precursor selection events }

\description{
  \code{precSelection} computes the number of selection events each
  precursor ions has undergone in an tandem MS experiment. This will be
  a function of amount of peptide loaded, chromatography efficiency,
  exclusion time,... and is useful when optimising and experimental
  setup. This function returns a named integer vector or length equal to
  the number of unique precursor MZ values in the original
  experiment. See \code{n} parameter to set the number of MZ significant
  decimals.

  \code{precSelectionTable} is a wrapper around \code{precSelection} and
  returns a table with the number of single, 2-fold, ... selection events.
}

\usage{
precSelection(object,n)
}

\arguments{
  \item{object}{ An instane of class \code{"\linkS4class{MSnExp}"}. }
  \item{n}{ The number of decimal places to round the precursor MZ
    to. Is passed to the \link{round} function. }
}


\value{
  A named integer in case of \code{precSelection} and a \code{table} for
  \code{precSelectionTable}. 
}


\author{
  Laurent Gatto <lg390@cam.ac.uk>
}

\examples{
precSelection(itraqdata)
precSelection(itraqdata,n=2)
precSelectionTable(itraqdata)
## only single selection event in this reduced exeriment
}
