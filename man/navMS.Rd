% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nav.R
\name{navMS}
\alias{navMS}
\alias{nextMS}
\alias{prevMS}
\title{Navigate an \code{MSnExp} object}
\usage{
navMS(i, object, msLevel, nav = c("nextMS", "prevMS"))

nextMS(...)

prevMS(...)
}
\arguments{
\item{i}{The name of index of the current spectrum}

\item{object}{The \code{MSnExp} object}

\item{msLevel}{The MS level of the next or previous spectrum. If
missing (default), then the level of the current spectrum is
used.}

\item{nav}{One of \code{"nextMS"} or \code{"prevMS"}, to obtain
the next or previous spectrum of level \code{msLevel}.}
}
\value{
An object of class \code{Spectrum1} or \code{Spectrum2},
    depending on the value of \code{msLevel} or \code{NULL}, of no
    spectrum is found.
}
\description{
Navigate an \code{MSnExp} object by moving to the next or previous
spectrum.
}
\examples{
f <- msdata::proteomics(full.names = TRUE)[1]
x <- readMSData2(f, centroided = c(FALSE, TRUE, FALSE))
sp <- "X009.1"
x[[sp]] ## curent MS3
nextMS(sp, x) ## next MS3
prevMS(sp, x) ## prev MS3
prevMS(sp, x, 2L) ## prev MS2
prevMS(sp, x, 1L) ## prev MS1
}
\author{
Laurent Gatto
}
