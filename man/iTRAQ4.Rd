\name{iTRAQ4}
\alias{iTRAQ4}
\alias{iTRAQ5}
\alias{iTRAQ8}
\alias{iTRAQ9}
\docType{data}

\title{ iTRAQ 4-plex set }

\description{
  This instance of class \code{"\linkS4class{ReporterIons}"} corresponds
  to the iTRAQ 4-plex set, i.e the 114, 115, 116 and 117 isobaric
  tags. In the iTRAQ5 data set, an unfragmented tag, i.e reporter and
  attached isobaric tag, is also included at MZ 145.
  These objects are used to plot the reporter ions of interest in an
  MSMS spectra (see \code{"\linkS4class{Spectrum2}"}) as well as for
  quantification (see \code{\link{quantify}}).
}

\usage{
data(iTRAQ4)
data(iTRAQ5)
data(iTRAQ8)
data(iTRAQ9)
}

\references{
  Ross PL, Huang YN, Marchese JN, Williamson B, Parker K, Hattan S,
  Khainovski N, Pillai S, Dey S, Daniels S, Purkayastha S, Juhasz P,
  Martin S, Bartlet-Jones M, He F, Jacobson A, Pappin DJ.
  "Multiplexed protein quantitation in Saccharomyces cerevisiae using
  amine-reactive isobaric tagging reagents."
  \emph{Mol Cell Proteomics}, 2004 Dec;3(12):1154-69.
  Epub 2004 Sep 22. PubMed PMID: 15385600.
}

\examples{
data(iTRAQ4)
iTRAQ4
data(iTRAQ5)
iTRAQ5
}

\seealso{
  \code{\link{TMT6}}.
}

\keyword{datasets}
