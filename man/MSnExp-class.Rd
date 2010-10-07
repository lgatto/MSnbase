\name{MSnExp-class}
\Rdversion{1.1}
\docType{class}
\alias{MSnExp-class}
\alias{class:MSnExp}
\alias{MSnExp}

\title{The "MSnExp" Class for MS Data And Meta-Data}

\description{
  The \code{MSnExp} class encapsulates data and meta-data for mass
  spectrometry experiments, as described in the \code{slots}
  section. Several data files (currently in \code{mzXML}) can be loaded
  together with the function \code{\link{readMzXMLData}}. 

  This class contains the virtual \code{"\linkS4class{eSet}"} class and
  adds MS-based proteomics specific attributes to describe the
  experiment, while making use of \code{"\linkS4class{eSet}"}'s existing
  slots for meta-data. 
}

\section{Slots}{
  \describe{
    \item{\code{spectra}:}{Object of class \code{"list"} of
      \code{"\linkS4class{Spectrum1}"} or \code{"\linkS4class{Spectrum2}"}
      objects - the data. }
    \item{\code{process}:}{Object of class
      \code{"\linkS4class{MSnProcess}"} that records all processing. }
    \item{\code{proteomicsData}:}{Object of class
      \code{"\linkS4class{MIAPE}"} that documents minimum information about
      the proteomics experiment. }
    \item{\code{description}:}{Object of class \code{"character"}
      describing the experiment in free form. }
    \item{\code{fromFile}:}{Object of class \code{"numeric"} of same
      length as \code{spectra} indicating the file of origin of each
      spectrum. }
    \item{\code{files}:}{Object of class \code{"character"} recording
      the name of the \code{mzXML} files that have been loaded by
      \code{\link{readMzXMLData}} to create the \code{MSnExp} instance. }
    \item{\code{featureData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing variables
      describing features (spectra in our case), e.g. identificaiton data,
      peptide sequence, identification score,... (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{featureData}} for
      more details. }
    \item{\code{phenoData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing 
      experimenter-supplied variables describing sample (i.e the
      individual tags for an labelled MS experiment) (inherited from 
      \code{"\linkS4class{eSet}"}). See \code{\link{phenoData}} for more
      details. }
    \item{\code{experimentData}:}{Object of class
      \code{"\linkS4class{MIAME}"}, containing details of experimental
      methods (inherited from \code{"\linkS4class{eSet}"}). See
      \code{\link{experimentData}} for more details.}
    \item{\code{protocolData}:}{Object of class
      \code{"\linkS4class{AnnotatedDataFrame}"} containing
      equipment-generated variables (inherited from
      \code{"\linkS4class{eSet}"}). See \code{\link{protocolData}} for
      more details. }
    \item{\code{assayData}:}{ not used here. }
    \item{\code{annotation}:}{ not used here. }
    \item{\code{.__classVersion__}:}{Object of class
      \code{"\linkS4class{Versions}"} describing the versions of R,
      \code{"linkS4class{Biobase}"}, \code{"\linkS4class{eSet}"} and
      \code{MSnExp} of the current instance. Intended for developer use
      and debugging (inherited from \code{"\linkS4class{eSet}"}). }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{eSet}"}, directly.
  Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
  Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}

%%\section{Methods}{
%%  No methods defined with class "MSnExp" in the signature.
%%}

\references{
  Information about the mzXML format as well converters from vendor
  specific formats to mzXML:
  \url{http://tools.proteomecenter.org/wiki/index.php?title=Formats:mzXML}. 
}
\author{
  Laurent Gatto <lg390@cam.ac.uk>
}
\seealso{
  \code{\link{readMzXMLData}} for loading \code{mzXML} data to generate
  an instance of \code{MSnExp}. See \code{"\linkS4class{eSet}"},
  \code{"\linkS4class{MSnProcess}"}, \code{"\linkS4class{MIAPE}"},
  \code{\link{featureData}}, \code{\link{phenoData}},
  \code{\link{experimentData}}, \code{"\linkS4class{MIAME}"}, 
  \code{\link{protocolData}}, for more details about the
  individual slots.
}
\examples{
mzxmlfile <- dir(system.file("extdata",package="MSnbase"),pattern="mzXML",full.names=TRUE)
msnexp <- readMzXMLData(mzxmlfile)
msnexp
}

\keyword{classes}
