\name{itraqdata}
\alias{itraqdata}
\docType{data}
\alias{itraqdata}
\alias{msnset}
\alias{msnset2}

\title{
  Example \code{MSnExp} and \code{MSnSet} data sets 
}

\description{

  \code{itraqdata} is and example data sets is an iTRAQ 4-plex
  experiment that has been run on an Orbitrap Velos instrument. It
  includes identification data in the feature data slot obtain from the
  Mascot search engine. It is a subset of an spike-in experiment where
  proteins have spiked in an \emph{Erwinia} background, as described in

  Karp et al. (2010), \emph{Addressing accuracy and precision issues in
    iTRAQ quantitation}, Mol Cell Proteomics. 2010 Sep;9(9):1885-97. Epub
  2010 Apr 10. (PMID 20382981).

  The spiked-in proteins in \code{itradata} are BSA and ENO and are
  present in relative abundances 1, 2.5, 5, 10 and 10, 5, 2.5, 1 in the
  114, 115, 116 and 117 reporter tags.

  The \code{msnset} object is produced by running the \code{quantify}
  method on the \code{itraqdata} experimental data, as detailed in the
  \code{\link{quantify}} example. This example data set is used in the
  MSnbase-demo vignette, available with \code{vignette("MSnbase-demo",
  package="MSnbase")}.

  The \code{msnset2} object is another example iTRAQ4 data that is used
  to demonstrate features of the package, in particular the \code{iPQF}
  feature aggregation method, described in \code{\link{iPQF}}. It
  corresponds to 11 proteins with spectra measurements from the original
  data set described by Breitwieser et al. (2011) \emph{General
  statistical modeling of data from protein relative expression isobaric
  tags}. J. Proteome Res., 10, 2758-2766.
  
}

\usage{itraqdata}

\examples{
data(itraqdata)
itraqdata

## created by
## msnset <- quantify(itraqdata, method = "trap", reporters = iTRAQ4)
data(msnset)
msnset

data(msnset2)
msnset2
}

\keyword{datasets}
