\name{TMT6}
\alias{TMT6}
\alias{TMT7}
\alias{TMT10}
\docType{data}

\title{ TMT 6-plex set }

\description{

  This instance of class \code{"\linkS4class{ReporterIons}"} corresponds
  to the TMT 6-plex set, i.e the 126, 127, 128, 129, 130 and 131
  isobaric tags. In the \code{TMT7} data set, an unfragmented tag, i.e
  reporter and attached isobaric tag, is also included at MZ 229. The
  \code{TMT10} instance corresponds to the 10-plex version.

  These objects are used to plot the reporter ions of interest in an
  MSMS spectra (see \code{"\linkS4class{Spectrum2}"}) as well as for
  quantification (see \code{\link{quantify}}).
  
}

\usage{
TMT6
TMT7
TMT10
}

\references{
  Thompson A, Sch\"{a}fer J, Kuhn K, Kienle S, Schwarz J, Schmidt G,
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
