\name{TMT6}
\alias{TMT6}
\alias{TMT6b}
\alias{TMT7}
\alias{TMT7b}
\alias{TMT10}
\alias{TMT10ETD}
\alias{TMT10HCD}
\alias{TMT11}
\alias{TMT11HCD}
\alias{TMT16}
\alias{TMT16HCD}
\docType{data}

\title{ TMT 6/10-plex sets }

\description{

  This instance of class \code{"\linkS4class{ReporterIons}"} corresponds
  to the TMT 6-plex set, i.e the 126, 127, 128, 129, 130 and 131
  isobaric tags. In the \code{TMT7} data set, an unfragmented tag, i.e
  reporter and attached isobaric tag, is also included at MZ 229. A
  second \code{TMT6b} has slightly different values.


  The \code{TMT10} instance corresponds to the 10-plex version. There
  are spectific HCD (\code{TMT10HCD}, same as \code{TMT10}) and ETD
  (\code{TMT10ETD}) sets.

  These objects are used to plot the reporter ions of interest in an
  MSMS spectra (see \code{"\linkS4class{Spectrum2}"}) as well as for
  quantification (see \code{\link{quantify}}).

}

\usage{
TMT6
TMT6b
TMT7
TMT7b
TMT10
TMT10ETD
TMT10HCD
TMT11
TMT11HCD
}

\references{
  Thompson A, Schäfer J, Kuhn K, Kienle S, Schwarz J, Schmidt G,
  Neumann T, Johnstone R, Mohammed AK, Hamon C.
  "Tandem mass tags: a novel quantification strategy for comparative
  analysis of complex protein mixtures by MS/MS."
  \emph{Anal Chem.} 2003 Apr 15;75(8):1895-904. \emph{Erratum} in:
  \emph{Anal Chem.} 2006 Jun 15;78(12):4235. Mohammed, A Karim A [added]
  and
  \emph{Anal Chem.} 2003 Sep 15;75(18):4942. Johnstone, R [added].
  PubMed PMID: 12713048.
}

\examples{
TMT6
TMT6[1:2]

TMT10

newReporter <- new("ReporterIons",
                   description="an example",
                   name="my reporter ions",
                   reporterNames=c("myrep1","myrep2"),
                   mz=c(121,122),
                   col=c("red","blue"),
                   width=0.05)
newReporter
}

\seealso{
  \code{\link{iTRAQ4}}.
}

\keyword{datasets}
