\name{itraqdata}
\alias{itraqdata}
\docType{data}

\title{
  Example \code{MSnExp} data set 
}

\description{
  This example data sets is an iTRAQ 4-plex experiment that
  has been run on an Orbitrap Velos instrument. It includes
  identification data in the feature data slot obtain from the Mascot
  search engine. 

  \code{itraqdata} is a subset of an spike-in experiment where proteins
  have spiked in an \emph{Erwinia} background, as described in
  Karp et al. (2010), \emph{Addressing accuracy and precision issues in
    iTRAQ quantitation}, Mol Cell Proteomics. 2010
  Sep;9(9):1885-97. Epub 2010 Apr 10. (PMID 20382981).
  The spiked-in proteins in \code{itradata} are BSA and ENO and are
  present in relative abundances 1, 2.5, 5, 10 and 10, 5, 2.5, 1 in the
  114, 115, 116 and 117 reporter tags.  
  
  This example data set is used in the MSnbase-demo vignette, available
  with \code{vignette("MSnbase-demo",package="MSnbase")}.
}

\usage{itraqdata}

\examples{
itraqdata
}

\keyword{datasets}
